# saveAll.py
from pymol import cmd

for obj in cmd.get_object_list():
    cmd.save(f"{obj}_align.pdb", obj)
