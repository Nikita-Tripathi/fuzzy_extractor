#!/bin/bash
while getopts 'c:e:' OPTION; do
    case "$OPTION" in
        c) 
            CODE=${OPTARG}
            # echo "Code is ${CODE}"
            ;;
        e)
            ERROR=${OPTARG}
            # echo "Error rate is ${ERROR}"
            ;;
        # ?)
        #     echo "script usage: $(basename \$0) -m path/to/message.src -s seed_integer -e error_rate_double"
        #     exit 1
        #     ;;
    esac
done

# Prepare 5G parity check matrix (612,198)
PCHK=./parity.pchk
if [ -f "$FILE" ]; then
    ./alist-to-pchk NR_1_4_18.alist parity.pchk
    ./make-gen parity.pchk gen.gen dense
fi

FILENAME="temp${RANDOM}.src"

echo ${CODE} >> ${FILENAME}

# Assuming input CODE is a path to a .rec file
# Decode the noisy message and store it to d.dec
./decode parity.pchk ${FILENAME} d.dec bsc "${ERROR}" prprp 250 2>&1

# Extract encoded message with errors and store this in r.rec
./extract gen.gen d.dec e.ext

cat e.ext

rm ${FILENAME}