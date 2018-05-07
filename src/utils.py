import os

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)
