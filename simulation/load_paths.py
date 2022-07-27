import os


def load_box_paths(user_path=None, parser_default='HPC'):

    if not user_path:
        user_path = os.path.expanduser('~')

    if 'mambrose' in user_path:
        home_path = os.path.join(user_path, 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'mAb_rtss_smc_comparison')
        datapath = os.path.join(home_path, 'data')
        projectpath = home_path
    else:
        home_path = os.path.join(user_path, 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'mAb_rtss_smc_comparison')
        datapath = os.path.join(home_path, 'data')
        projectpath = home_path

    return datapath, projectpath
