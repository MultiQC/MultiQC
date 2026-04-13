"""MultiQC Utility functions, used in a variety of places."""

import array
import json
import logging
import math
import shutil
import sys
import time
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
from pydantic import BaseModel

logger = logging.getLogger(__name__)


def rmtree_with_retries(
    path: Union[str, Path, None],
    _logger: Optional[logging.Logger] = None,
    max_retries: int = 10,
):
    """
    Robustly tries to delete paths.
    Retries several times (with increasing delays) if an OSError
    occurs.  If the final attempt fails, the Exception is propagated
    to the caller.
    """
    if path is None or not Path(path).exists():
        return

    for i in range(max_retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            if _logger:
                _logger.info(f"Unable to remove path: {path}")
                _logger.info(f"Retrying after {i**2} seconds")
            else:
                print(f"Unable to remove path: {path}", file=sys.stderr)
                print(f"Retrying after {i**2} seconds", file=sys.stderr)
            time.sleep(i**2)

    # Final attempt, pass any Exceptions up to caller.
    shutil.rmtree(path)


def strtobool(val: Any) -> bool:
    """
    Replaces deprecated https://docs.python.org/3.9/distutils/apiref.html#distutils.util.strtobool
    The deprecation recommendation is to re-implement the function https://peps.python.org/pep-0632/

    ------------------------------------------------------------

    Convert a string representation of truth to true (1) or false (0).

    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val_str = str(val).lower()
    if val_str in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val_str in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")


def replace_defaultdicts(data: Any) -> Any:
    """
    Recursively replace dict-likes as dicts for nice yaml representation.
    """

    def _replace(obj: Any) -> Any:
        if isinstance(obj, (defaultdict, OrderedDict, dict)):
            return {k: _replace(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_replace(v) for v in obj]
        elif isinstance(obj, set):
            return {_replace(v) for v in obj}
        elif isinstance(obj, tuple):
            return tuple(_replace(v) for v in obj)
        return obj

    return _replace(data)


def dump_json(data, filehandle=None, **kwargs):
    """
    Recursively replace non-JSON-conforming NaNs and lambdas with None.
    Note that a custom JSONEncoder would not work for NaNs:
    https://stackoverflow.com/a/28640141
    """

    def replace_nan(obj):
        """
        Recursively replace NaNs and Infinities with None
        """
        # Do checking in order of likelihood of occurrence
        if isinstance(obj, float):
            if math.isnan(obj) or math.isinf(obj):
                return None
            return obj
        if isinstance(obj, (tuple, set)):
            # JSON only knows list so convert tuples and sets to list.
            obj = list(obj)
        if isinstance(obj, list):
            for i, item in enumerate(obj):
                if isinstance(item, float) and (math.isnan(item) or math.isinf(item)):
                    obj[i] = None
                elif isinstance(item, (dict, list, tuple, set)):
                    obj[i] = replace_nan(item)
            return obj
        if isinstance(obj, dict):
            for key, value in obj.items():
                if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
                    obj[key] = None
                elif isinstance(value, (dict, list, tuple, set)):
                    obj[key] = replace_nan(value)
            return obj
        return obj

    class JsonEncoderWithArraySupport(json.JSONEncoder):
        """
        Encode array.array instances to list. Use the default method
        for this as it gets called only when an array instance is encountered
        and is then immediately serialized into a string. This saves memory
        compared to unpacking all arrays to list at once.
        """

        def default(self, o):
            if isinstance(o, array.array):
                return replace_nan(o.tolist())
            if callable(o):
                return None
            if isinstance(o, BaseModel):  # special handling for pydantic models
                return o.model_dump_json()
            return super().default(o)

    if filehandle:
        try:
            json.dump(replace_nan(data), filehandle, cls=JsonEncoderWithArraySupport, **kwargs)
        except UnicodeEncodeError:
            # If we get an error, try to write the JSON string directly
            # First convert to a string with proper encoding
            json_str = json.dumps(replace_nan(data), cls=JsonEncoderWithArraySupport, ensure_ascii=False, **kwargs)
            # Write the string directly to avoid encoding issues
            if hasattr(filehandle, "write"):
                filehandle.write(json_str)
            else:
                with open(filehandle, "w", encoding="utf-8") as f:
                    f.write(json_str)
    else:
        return json.dumps(replace_nan(data), cls=JsonEncoderWithArraySupport, **kwargs)


def is_running_in_notebook() -> bool:
    try:
        from IPython import get_ipython  # type: ignore

        if "IPKernelApp" in get_ipython().config:
            return True
    except (ImportError, AttributeError):
        pass
    return False


def compress_number_lists_for_json(obj):
    """
    Take an object that should be JSON and compress all the lists of integer
    and lists of float as array.array. This saves space and the arrays can
    easily be converted back again, using the dump_json function above.

    The technical explanation:
    A python list is an array of pointers to python objects:
    {
        list metadata including length
        a pointer to an array of pointers: [
            PyObject *
            PyObject *
            etc.
        ]
    }
    A python float is very simple and takes 24 bytes.
    {
        PyTypeObject *type
        Py_ssize_t refcount
        double the actual floating point.
    }
    A python integer is slightly more complicated, but similar to the float. It
    takes 28 bytes for 32-bit data, 32 bytes for 64-bit, 36 bytes for 96-bit etc.

    An array.array is more simple.
    {
        array metadata including length
        a pointer to an array of machine values: [
            double,
            double,
            double,
            etc.
        ]
        more metadata
    }
    Using 8-byte machine values rather than Python objects saves thus
    24 bytes per float.
    """
    if isinstance(obj, (list, tuple)):
        try:
            # Try integer list first, because it does not accept floats.
            return array.array("q", obj)
        except TypeError:
            pass
        try:
            return array.array("d", obj)
        except TypeError:
            return [compress_number_lists_for_json(v) for v in obj]
    if isinstance(obj, dict):
        return {k: compress_number_lists_for_json(v) for k, v in obj.items()}
    return obj


def update_dict(
    target: Dict[Any, Any],
    source: Dict[Any, Any],
    none_only: bool = False,
    add_in_the_beginning: bool = False,
):
    """
    Recursively updates nested dict d from nested dict u

    >>> update_dict({"cutadapt": {"fn": "old", "fn2": "old2"}}, {"cutadapt": {"fn": "new"}})
    {'cutadapt': {'fn': 'new', 'fn2': 'old2'}}
    >>> update_dict({"cutadapt": [{"fn": "old"}]}, {"cutadapt": {"fn": "new"}})
    {'cutadapt': {'fn': 'new'}}
    >>> update_dict({"existing": "v1"}, {"new": "v2"})
    {'existing': 'v1', 'new': 'v2'}
    >>> update_dict({"existing": "v1"}, {"new": "v2"}, add_in_the_beginning=True)
    {'new': 'v2', 'existing': 'v1'}
    """
    for key, src_val in source.items():
        if isinstance(src_val, dict) and key in target and isinstance(target[key], dict):
            target[key] = update_dict(target[key], src_val, none_only=none_only)
        else:
            if not none_only or target.get(key) is None:
                if isinstance(src_val, list):
                    target[key] = src_val.copy()
                else:
                    if add_in_the_beginning:
                        target = {key: src_val, **target}
                    else:
                        target[key] = src_val
    return target


def scipy_pdist(X: np.ndarray) -> np.ndarray:
    """Calculate pairwise (euclidean) distances between observations in X.

    Reimplements scipy.spatial.distance.pdist to avoid heavy scipy dependency.

    Args:
        X: Array of shape (m,n) containing m observations in n dimensions

    Returns:
        Array of shape ((m * (m-1)) // 2,) containing condensed distance matrix
    """
    m, n = X.shape
    # Initialize output array of correct size for condensed distance matrix
    out = np.zeros((m * (m - 1)) // 2)
    k = 0

    # Calculate pairwise distances
    for i in range(m - 1):
        for j in range(i + 1, m):
            out[k] = np.sqrt(np.sum((X[i] - X[j]) ** 2))
            k += 1

    return out


def scipy_hierarchy_linkage(distances: np.ndarray, method: str = "complete") -> np.ndarray:
    """Perform hierarchical clustering using the specified linkage method.

    Reimplements scipy.hierarchy.linkage to avoid heavy scipy dependency.

    Args:
        distances: Condensed distance matrix from pdist
        method: Linkage method ('single', 'complete', 'average', 'weighted')

    Returns:
        Array of shape (n-1, 4) representing the linkage matrix where n is the number of original observations.
        Each row has format [idx1, idx2, distance, cluster_size]
    """
    if method not in ["single", "complete", "average", "weighted"]:
        raise ValueError(f"Unsupported linkage method: {method}")

    # Convert condensed distance matrix to full matrix for easier manipulation
    n = int((1 + np.sqrt(1 + 8 * len(distances))) / 2)
    dist_matrix = np.zeros((n, n))
    idx = np.triu_indices(n, k=1)
    dist_matrix[idx] = distances
    dist_matrix = dist_matrix + dist_matrix.T  # type: ignore

    # Initialize arrays for clustering
    n_clusters = n
    active_nodes = list(range(n_clusters))
    cluster_sizes = np.ones(n_clusters, dtype=int)
    linkage_matrix = np.zeros((n_clusters - 1, 4))

    for i in range(n_clusters - 1):
        # Find minimum distance between clusters
        valid_distances = dist_matrix[np.ix_(active_nodes, active_nodes)]
        np.fill_diagonal(valid_distances, np.inf)
        min_idx = np.unravel_index(np.argmin(valid_distances), valid_distances.shape)
        cluster1, cluster2 = active_nodes[min_idx[0]], active_nodes[min_idx[1]]

        # Record merge in linkage matrix
        linkage_matrix[i] = [
            min(cluster1, cluster2),
            max(cluster1, cluster2),
            dist_matrix[cluster1, cluster2],
            cluster_sizes[cluster1] + cluster_sizes[cluster2],
        ]

        # Create new cluster
        new_cluster_idx = n_clusters + i
        cluster_sizes = np.append(cluster_sizes, cluster_sizes[cluster1] + cluster_sizes[cluster2])  # type: ignore

        # Calculate distances to new cluster based on chosen method
        remaining_clusters = [x for x in active_nodes if x not in (cluster1, cluster2)]
        new_distances = np.zeros(len(remaining_clusters))

        for j, cluster in enumerate(remaining_clusters):
            dist1 = dist_matrix[cluster1, cluster]
            dist2 = dist_matrix[cluster2, cluster]

            if method == "single":
                new_dist = min(dist1, dist2)
            elif method == "complete":
                new_dist = max(dist1, dist2)
            elif method == "average":
                new_dist = (dist1 * cluster_sizes[cluster1] + dist2 * cluster_sizes[cluster2]) / (
                    cluster_sizes[cluster1] + cluster_sizes[cluster2]
                )
            else:  # weighted
                new_dist = (dist1 + dist2) / 2

            new_distances[j] = new_dist

        # Update distance matrix
        dist_matrix = np.vstack((dist_matrix, np.zeros(dist_matrix.shape[1])))  # type: ignore
        dist_matrix = np.hstack((dist_matrix, np.zeros((dist_matrix.shape[0], 1))))  # type: ignore
        dist_matrix[new_cluster_idx, remaining_clusters] = new_distances  # type: ignore
        dist_matrix[remaining_clusters, new_cluster_idx] = new_distances  # type: ignore

        # Update active nodes
        active_nodes.remove(cluster1)
        active_nodes.remove(cluster2)
        active_nodes.append(new_cluster_idx)

    return linkage_matrix


def scipy_hierarchy_leaves_list(Z: np.ndarray) -> List[int]:
    """Return the leaf nodes in the order they appear in the dendrogram.

    Reimplements scipy.hierarchy.leaves_list to avoid heavy scipy dependency.

    Args:
        Z: The linkage matrix from scipy_hierarchy_linkage

    Returns:
        List of original observation indices in the order they appear in the dendrogram
    """
    n = int(Z.shape[0] + 1)
    # Initialize list of current clusters with original observations
    clusters = [[i] for i in range(n)]

    # Process merges in order
    for i in range(len(Z)):
        # Get the clusters being merged
        cluster1 = int(Z[i, 0])
        cluster2 = int(Z[i, 1])

        # Merge clusters
        new_cluster = clusters[cluster1] + clusters[cluster2]
        clusters.append(new_cluster)

    # Return the final ordering (last cluster created)
    return clusters[-1]
