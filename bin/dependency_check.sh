
ret=0

# check to see if the appropriate programs are found...
# picard
PicardCommandLine -h 2>/dev/null
if [ $? -ne 1 ]; then # don't ask me why it returns 1.
    echo "Please install Picard!"
    ret=1
fi

# GNU's parallel
parallel -h >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install GNU parallel!"
    ret=1
fi

# GATK (4.*)
gatk --help >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install GATK!"
    ret=1
fi

#what's hap! (physical phasing)
whatshap --help >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install whatshap!!"
    ret=1
fi


bcftools -h >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install bcftools!"
    ret=1
fi

tabix -h 2>/dev/null
if [ $? -ne 1 ]; then # check for the correct return value (1 in this case... 'cause they're annoying)
    echo "Please install tabix!!"
    ret=1
fi

shapeit4 &> /dev/null
if [ $? -ne 1 ]; then # check for the return val
    echo "Please install shapeit4!!"
    ret=1
fi

exit $ret


