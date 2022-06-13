from setuptools import setup
import sys

def print_success(message=None):
		if message is None:
			sys.stdout.write("OK\n")
		else:
			message = str(message)
			sys.stdout.write(message)
			if not message.endswith("\n"):
				sys.stdout.write("\n")
		sys.stdout.flush()
if __name__ == "__main__":
	try:
		setup()
		print_success('While KiMoPack is a python library, we facilitate its use with Jupyter notebooks. For the typical analysis tasks we have developed a series of Notebooks that guide through the tasks.\n These notebooks can be found in the installation folder under "Workflow_tools" or can be downloaded from https://github.com/erdzeichen/KiMoPack/tree/main/Workflow_tools. Please copy one of these notebooks into your data analysis folder and rename them to create a analysis log of your session. For more information please see the publication https://doi.org/10.1021/acs.jpca.2c00907, the tutorial videos, or the tutorial notebooks under https://github.com/erdzeichen/KiMoPack/tree/main/Tutorial_Notebooks_for_local_use. ')
	except Exception as e:
		print_success(message)