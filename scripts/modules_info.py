"""
Loop over the modules code and print summary about each module:
 - how many plots are used by each module, split by the plot type.

Potentially will print more per-module info, as well as aggregate stats.
"""

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict


def find_plot_calls(module_path: Path) -> Dict[str, int]:
    call_by_plot_type: Dict[str, int] = defaultdict(int)
    with module_path.open() as f:
        for line in f:
            if ".plot(" in line:
                if m := re.search(r"(\w+)\.plot\(", line):
                    call_by_plot_type[m.group(1)] += 1
            if "general_stats_addcols(" in line:
                call_by_plot_type["general_stats_table"] = 1
    return call_by_plot_type


def main():
    modules_dir = Path(__file__).parent.parent / "multiqc/modules"
    stat_by_plot_type_by_module: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for mod_dir in modules_dir.iterdir():
        if mod_dir.is_dir():
            mod_name = mod_dir.name
            for mod_file in mod_dir.iterdir():
                if mod_file.suffix == ".py":
                    for plot_type, count in find_plot_calls(mod_file).items():
                        stat_by_plot_type_by_module[mod_name][plot_type] += count

    sorted_stats = dict(sorted(stat_by_plot_type_by_module.items(), key=lambda x: sum(x[1].values()), reverse=True))
    for mod_name, stats_by_plot_type in sorted_stats.items():
        print(mod_name)
        for plot_type, count in stats_by_plot_type.items():
            print(f"  {plot_type}: {count}")


if __name__ == "__main__":
    main()
