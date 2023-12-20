# Disassembly

The [idea](/paper/idea.md).

TODO:

1. figure out the filtering problem with proteolysis simulation
2. make the weight-updating better. 
    - Take exoprotease into account. 
    - Set higher probability to ABCD - ABC if ABC and D is present. 
    - Increase weight ABCD - ABC more if BOTH ABCD and ABC is of high abundance. Right now it only takes abundance of ABC into account.