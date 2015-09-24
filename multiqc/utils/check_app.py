import os

def check_app():
    """Check is OSX app is running
    
        Return
            osx_app (bool): If app is running or not
    """
    # Are we running as an OSX App?
    if os.environ.get('MULTQC_IS_APP') is not None:
        osx_app = True
        logging.basicConfig(stream=sys.stdout, format='<li>%(message)s</li>', level=20)
        with open(os.path.join(os.getcwd(), 'multiqc_app_header.html'), 'r') as f:
            print(f.read())
    else:
        osx_app = False
    
    return osx_app
