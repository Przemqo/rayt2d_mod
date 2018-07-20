 #!/bin/bash
 gcc -I./include/ -L./lib/ \
 ./source/rayt2d.c -lpar -lcwp -lm \
 -o /home/przemek/Desktop/for_przemyslaw/compiled/rayt2d2

exit 0
#Remove path from last line, leave rayt2d_mod
