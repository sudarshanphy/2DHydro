import os
import sys
import shutil

def parse_arguments():
    """
    Parse command-line arguments in the form +key=value.
    Returns a dictionary of arguments with default values.

    Defaults:
        make: 'linux'
        objdir: 'build'
        runtype: 'hd'
        offload: 'cpu'
    """
    # Define default values for all possible arguments
    defaults = {
        'make': 'linux',
        'objdir': 'build',
        'runtype': 'hd',
        'offload': 'cpu'
    }

    # Parse command-line arguments, overriding defaults
    for arg in sys.argv[1:]:
        if arg.startswith('+') and '=' in arg:
            key, value = arg[1:].split('=', 1)
            if key in defaults:  # Only accept known keys
                defaults[key] = value
    return defaults

def validate_paths(args):
    """
    Validate that required files and directories exist.
    Returns a tuple of validated paths and module implementations.
    """
    FilePath = "./"
    # Validate Makefile location
    site_dir = os.path.join(FilePath+'sites', args['make'])
    site_makefile = os.path.join(site_dir, 'Makefile')
    if not os.path.exists(site_makefile):
        raise FileNotFoundError(f"Makefile not found in site directory: {site_dir}")
    
    # Create object directory if needed
    objdir = FilePath+args['objdir']
    if os.path.exists(objdir):
        shutil.rmtree(objdir)
    os.makedirs(objdir)

    return site_makefile, objdir


def copy_files(objdir, args):
    src_dir = "src"
    dest_dir = objdir
    for root, dirs, files in os.walk(src_dir):
        # Create corresponding destination directory
        rel_path = os.path.relpath(root, src_dir)
        dest_path = os.path.join(objdir)
        #print(root)
        #print(dirs)
        #print(files)
        #print(rel_path)
        #print(dest_path)
        #print("------------------")
        if (root == src_dir):
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path, file)
                #print(src_file, dest_file)
                shutil.copy(src_file, dest_path)
        if ((args['runtype'] == "hd") and (rel_path == "HD_h")):
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path,file)
                shutil.copy(src_file, dest_file)
        if ((args['runtype'] == "mhd") and (rel_path == "MHD_h")):
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path,file)
                shutil.copy(src_file, dest_file)
        if ((args['offload'] == "gpu") and (rel_path == "GPU")):
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path,file)
                shutil.copy(src_file, dest_file)
           
    return

       

def main():
    try:
        # Parse arguments with defaults
        args = parse_arguments()
        print("Configuration settings:")
        for k, v in args.items():
            print(f"  {k}: {v}")

        # Validate paths and get module implementations
        site_makefile, objdir = validate_paths(args)

        # Copy Makefile
        shutil.copy(site_makefile, os.path.join(objdir, 'Makefile'))

        copy_files(objdir, args)

        print("\nBuild configuration completed successfully!")
        print(f"Object directory: {objdir}")
        
    except Exception as e:
        print(f"\nError: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

