import os
import shutil
import tempfile

project_temp = os.environ['PROJECT_TEMP']

def temp_dir(name, args):
    '''
    Returns a temporary directory that will be deletd
    when python finishes.
    If used in a context manager, then this is just the string
    path of that directory and will be cleaned up when exiting
    the context.
    Otherwise will be a TemporaryDirectory object that has
    the attribution .name and the method .cleanup() to
    delete the directory.
    '''
    dname = name + '-'
    first = True
    for k, v in vars(args).items():
        if not first:
            dname += ','
        first = False
        dname += k.replace('/', '_') + '=' + str(v).replace('/', '_')
    dname += '-'
    return tempfile.TemporaryDirectory(prefix=dname, dir=project_temp)

def move_files(temp_dir, permanent_loc):
    assert 'lustre' in temp_dir
    for fname in os.listdir(temp_dir):
        shutil.move(os.path.join(temp_dir, fname), permanent_loc)

