import os
import shutil
import tempfile

ukb = os.environ['UKB']
project_temp = ukb + '/scratch'

def temp_dir(outdir):
    '''
    Returns a temporary directory that will be deleted
    when python finishes.
    If used in a context manager, then this is just the string
    path of that directory and will be cleaned up when exiting
    the context.
    Otherwise will be a TemporaryDirectory object that has
    the attribution .name and the method .cleanup() to
    delete the directory.
    '''
    if outdir.startswith(ukb):
        outdir = outdir[len(ukb):]
        if outdir.startswith('/'):
            outdir = outdir[1:]
    split = outdir.split('/')
    if len(split) > 1:
        pre_dir = project_temp + '/' + '/'.join(split[:-1])
        os.makedirs(pre_dir, exist_ok=True)
    else:
        pre_dir = project_temp
   
    print(split[-1], pre_dir)
    tdir = tempfile.TemporaryDirectory(prefix=split[-1], dir=pre_dir)
    print(f"Working with tempdir {tdir.name}")
    return tdir

def move_files(temp_dir, permanent_loc):
    for fname in os.listdir(temp_dir):
        shutil.move(os.path.join(temp_dir, fname), permanent_loc)

