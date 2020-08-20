rm -rf build dist *.egg-info
conda install setuptools wheel twine -y
conda update setuptools wheel twine -y

python setup.py sdist bdist_wheel
twine upload dist/* --verbose