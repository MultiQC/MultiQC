import logging

def get_sorted_modules(multiqc_submodules, module_order):
    """Check if all modules have a order
    
        There is an option to display the modules in a certain order in the 
        report. For example it couls be visually appealing to display qc data
        in the same order as the data have been processed.
        
        The order is specified in a list called module_order.
        This function checks if all available modules have a order.
        If not it logs a warning.
        
        The module will still be used but are just added to the end of the list.
        
        Args:
            multiqc_submodules (iterator): These are all available modules
            module_order (list): List with the order of the modules
        
        Returns:
            sorted_modules : A list with all modules
    """
    logger = logging.getLogger(__name__)
    
    logger.debug("Adding modules in sorted order to sorted_modules")
    sorted_modules = [module for module in module_order if module in multiqc_submodules]
    
    if len(multiqc_submodules) > len(sorted_modules):
        logger.warning("Modules missing from order declaration: {}".format(
            ', '.join([m for m in multiqc_submodules if m not in module_order])))
    
        logger.debug("Adding unordered modules to sorted_modules")
        sorted_modules.extend([module for module in multiqc_submodules if module_order not in module_order])
    
    return sorted_modules
    