 #!/bin/bash
 gcc -I./include/ -L./lib/ \
 ./source/rayt2d_mod.c -lpar -lcwp -lm \
 -o /home/przemek/Desktop/for_przemyslaw/compiled/rayt2d_mod

exit 0
#Remove path from last line, leave rayt2d_mod
