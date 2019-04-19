import os
for higherZ in range(70,84):
    filename = 'sh_iso_'+str(higherZ)+'py'
    print(filename)
    writestring ="""#!/bin/sh
set -x

#    Get isodata
$GRASP/bin/iso <<EOF
""" + str(higherZ) + """
171
n
170.9363258
0.5
0.49367
1.6
EOF

cat isodata"""

    f = open(filename,'w')
    f.write(writestring)
    f.close()
