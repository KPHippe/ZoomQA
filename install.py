import pathlib

path_file = open("./script/paths.py", 'r').readlines()
install_path = pathlib.Path(__file__).parent.resolve()

# use the python script to change the python script
for i in range(len(path_file)):
    # print(path_file[i])
    if "self.sw_install =" in path_file[i]:
        path_file[i] = f"        self.sw_install = '{install_path}'\n"


with open("./script/paths.py", 'w+') as f:
    f.write("".join(path_file))
