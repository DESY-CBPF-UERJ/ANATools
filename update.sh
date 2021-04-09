git fetch
git merge origin/main
python3 setup.py sdist
pip3 install ./dist/anatools-0.1.0.tar.gz