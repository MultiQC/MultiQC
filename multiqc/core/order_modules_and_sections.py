import logging

from multiqc import config, report
from multiqc.core.file_search import include_or_exclude_modules

logger = logging.getLogger(__name__)


def order_modules_and_sections():
    """
    Finalise modules and sections: place in the write order, add special-case modules.
    """
    # Importing here to avoid circular imports
    from multiqc.modules.software_versions import MultiqcModule as SoftwareVersionsModule
    from multiqc.modules.profile_runtime import MultiqcModule as ProfileRuntimeModule

    # First, remove the special-case modules that we want to re-add at the end, in case if they were
    # already added by a previous call of multiqc.write_report in an interactive session.
    report.modules = [
        m
        for m in report.modules
        if (not isinstance(m, SoftwareVersionsModule) and not isinstance(m, ProfileRuntimeModule))
    ]

    # In case if user passed exclude_modules or include_modules again:
    mod_ids = include_or_exclude_modules([mod.id for mod in report.modules])
    for mod in report.modules:
        if mod.id not in mod_ids:
            mod.hidden = True

    # Add section for software versions if any are found
    if not config.skip_versions_section and report.software_versions:
        report.modules.append(SoftwareVersionsModule())

    # Special-case module if we want to profile the MultiQC running time
    if config.profile_runtime:
        report.modules.append(ProfileRuntimeModule())

    # Sort the report module output if we have a config
    if len(config.report_section_order) > 0:
        section_id_order = {}
        idx = 10
        for mod in reversed(report.modules):
            section_id_order[mod.anchor] = idx
            idx += 10
        for anchor, ss in config.report_section_order.items():
            if anchor not in section_id_order.keys():
                logger.debug(f"Reordering sections: anchor '{anchor}' not found.")
                continue
            if ss.get("order") is not None:
                section_id_order[anchor] = ss["order"]
            if ss.get("after") in section_id_order.keys():
                section_id_order[anchor] = section_id_order[ss["after"]] + 1
            if ss.get("before") in section_id_order.keys():
                section_id_order[anchor] = section_id_order[ss["before"]] - 1
        sorted_ids = sorted(section_id_order.keys(), key=lambda k: section_id_order[k])
        report.modules = [mod for i in reversed(sorted_ids) for mod in report.modules if mod.anchor == i]

    # Sort the report sections if we have a config
    # Basically the same as above, but sections within a module
    if len(config.report_section_order) > 0:
        # Go through each module
        for midx, mod in enumerate(report.modules):
            section_id_order = dict()
            # Get a list of the section anchors
            idx = 10
            for s in mod.sections:
                section_id_order[s.anchor] = idx
                idx += 10
            # Go through each section to be reordered
            for anchor, ss in config.report_section_order.items():
                # Section to be moved is not in this module
                if anchor not in section_id_order.keys():
                    logger.debug(f"Reordering sections: anchor '{anchor}' not found for module '{mod.name}'.")
                    continue
                if ss == "remove":
                    section_id_order[anchor] = False
                    continue
                if ss.get("order") is not None:
                    section_id_order[anchor] = ss["order"]
                if ss.get("after") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["after"]] + 1
                if ss.get("before") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["before"]] - 1
            # Remove module sections
            section_id_order = {s: o for s, o in section_id_order.items() if o is not False}
            # Sort the module sections
            sorted_ids = sorted(section_id_order.keys(), key=lambda k: section_id_order[k])
            report.modules[midx].sections = [s for i in sorted_ids for s in mod.sections if s.anchor == i]
