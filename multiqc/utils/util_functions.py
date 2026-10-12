"""MultiQC Utility functions, used in a variety of places."""

import array
import json
import logging
import math
import shutil
import sys
import time
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
        if isinstance(obj, dict):
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

        ipython = get_ipython()
        if ipython is not None and "IPKernelApp" in ipython.config:
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
    # Contiguous rows make the summation order, and so tie-breaking, independent of the input memory layout
    X = np.ascontiguousarray(X)
    m, n = X.shape
    # Initialize output array of correct size for condensed distance matrix
    out = np.zeros((m * (m - 1)) // 2)
    k = 0

    # Calculate pairwise distances
    for i in range(m - 1):
        out[k : k + m - i - 1] = np.sqrt(np.sum((X[i] - X[i + 1 :]) ** 2, axis=1))
        k += m - i - 1

    return out


def scipy_hierarchy_linkage(distances: np.ndarray, method: str = "complete") -> np.ndarray:
    """Perform hierarchical clustering using the specified linkage method.

    Reimplements scipy.hierarchy.linkage to avoid heavy scipy dependency.

    Each step merges the closest pair of clusters, taking the smallest cluster ids on ties. A merged
    cluster reuses the matrix slot of its smaller id, and every cluster caches its closest partner
    among clusters with a larger id. When that partner is merged away, the cached distance stays as
    a lower bound and the row is rescanned only if it comes up as the closest.

    Args:
        distances: Condensed distance matrix from pdist
        method: Linkage method ('single', 'complete', 'average', 'weighted')

    Returns:
        Array of shape (n-1, 4) representing the linkage matrix where n is the number of original observations.
        Each row has format [idx1, idx2, distance, cluster_size]
    """
    if method not in ["single", "complete", "average", "weighted"]:
        raise ValueError(f"Unsupported linkage method: {method}")

    n = int((1 + np.sqrt(1 + 8 * len(distances))) / 2)
    if len(distances) != n * (n - 1) // 2:
        raise ValueError(f"Expected a condensed distance matrix, got {len(distances)} distances")
    dist_matrix = np.full((n, n), np.inf)
    start = 0
    for slot in range(n - 1):
        stop = start + n - slot - 1
        dist_matrix[slot, slot + 1 :] = distances[start:stop]
        dist_matrix[slot + 1 :, slot] = distances[start:stop]
        start = stop

    slot_ids = np.arange(n)
    slot_sizes = np.ones(n, dtype=int)
    # Active slots, ascending by cluster id
    order = np.arange(n)
    nearest = np.full(n, -1)
    nearest_dist = np.full(n, np.inf)
    stale = np.zeros(n, dtype=bool)
    for slot in range(n - 1):
        partner = slot + 1 + int(np.argmin(dist_matrix[slot, slot + 1 :]))
        nearest[slot] = partner
        nearest_dist[slot] = dist_matrix[slot, partner]

    linkage_matrix = np.zeros((n - 1, 4))
    for i in range(n - 1):
        # argmin returns the first NaN or the first minimum, so ties go to the smallest cluster id
        pos = int(np.argmin(nearest_dist[order]))
        while stale[order[pos]]:
            slot = order[pos]
            row = dist_matrix[slot, order[pos + 1 :]]
            best = int(np.argmin(row))
            nearest[slot] = order[pos + 1 + best]
            nearest_dist[slot] = row[best]
            stale[slot] = False
            pos = int(np.argmin(nearest_dist[order]))
        slot1 = order[pos]
        if nearest_dist[slot1] == np.inf:
            raise ValueError("Cannot cluster observations at infinite distance")
        slot2 = nearest[slot1]
        size1, size2 = slot_sizes[slot1], slot_sizes[slot2]
        linkage_matrix[i] = [slot_ids[slot1], slot_ids[slot2], dist_matrix[slot1, slot2], size1 + size2]

        rest = order[(order != slot1) & (order != slot2)]
        dist1 = dist_matrix[slot1, rest]
        dist2 = dist_matrix[slot2, rest]
        # Unlike np.minimum and np.maximum, np.where keeps dist1 when dist2 is NaN
        if method == "single":
            new_distances = np.where(dist2 < dist1, dist2, dist1)
        elif method == "complete":
            new_distances = np.where(dist2 > dist1, dist2, dist1)
        elif method == "average":
            new_distances = (dist1 * size1 + dist2 * size2) / (size1 + size2)
        else:  # weighted
            new_distances = (dist1 + dist2) / 2

        dist_matrix[slot1, rest] = new_distances
        dist_matrix[rest, slot1] = new_distances
        slot_ids[slot1] = n + i
        slot_sizes[slot1] = size1 + size2
        order = np.append(rest, slot1)

        current = nearest_dist[rest]
        closer = (new_distances < current) | (np.isnan(new_distances) & ~np.isnan(current))
        stale[rest[(nearest[rest] == slot1) | (nearest[rest] == slot2)]] = True
        stale[rest[closer]] = False
        nearest[rest[closer]] = slot1
        nearest_dist[rest[closer]] = new_distances[closer]
        nearest[slot1] = -1
        nearest_dist[slot1] = np.inf

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
