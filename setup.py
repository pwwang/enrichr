from setuptools import setup, find_packages

# get version
from os import path
verfile = path.join(path.dirname(__file__), 'enrichr.py')
with open(verfile) as vf:
    VERSION = vf.readline().split('=')[1].strip()[1:-1]

setup (
	name             = 'enrichr',
	version          = VERSION,
	description      = "Python wrapper for Enrichr APIs.",
	url              = "https://github.com/pwwang/enrichr",
	author           = "pwwang",
	author_email     = "pwwang@pwwang.com",
	license          = "MIT License",
	packages         = find_packages(),
	install_requires = [
		'matplotlib', 'requests'
    ],
)
