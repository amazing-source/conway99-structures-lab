import json, subprocess, sys
cfg = open(sys.argv[1]).read().strip()
sys.argv = ['fullblock_rigid.py', cfg] + sys.argv[2:]
exec(open('fullblock_rigid.py').read())
