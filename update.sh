git fetch
git merge origin/main
python3 setup.py sdist
pip3 install ./dist/anatools-0.2.0.tar.gz