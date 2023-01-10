#!/bin/bash
while getopts 's:e:' OPTION; do
    case "$OPTION" in
        s)
            SEED=${OPTARG}
            # echo "Seed is ${SEED}"
            ;;
        e)
            ERROR=${OPTARG}
            # echo "Error rate is ${ERROR}"
            ;;

    esac
done

# Prepare 5G parity check matrix (612,198)
PCHK=./parity.pchk
if [ -f "$FILE" ]; then
    ./alist-to-pchk NR_1_4_18.alist parity.pchk
    ./make-gen parity.pchk gen.gen dense
fi

# FILENAME="temp${RANDOM}.src"

# Put message into a temp file
# echo ${MESSAGE} >> ${FILENAME}

# Encode the message and send it to transmit
# Assuming SEED is an integer and ERROR is a double
# Populate encoded message with errors and send to stdout
./transmit e.enc r.rec "${SEED}" bsc "${ERROR}"

# Delete temp file
# rm ${MESSAGE}