# IMPORTS
from math import ceil
import numpy as np
from secrets import randbits
from subprocess import check_output
import galois
from PIL import Image
import multiprocessing
import hashlib

import sys, time, random
import mx_par



# RUNS ON PYTHON 3.8.9 +
#   (on my system its in /usr/bin/python3, with /usr/bin/pip3 package manager)

class FuzzyExtractor:

    def __init__(self, k=43, l=1000000, err=0.11, ecc_len=1224, lbd=128, xi = 128, t=12, file_prefix=''):
        '''
        Specify subsample size, number of lockers. and LPN encryption error rate
        Default parameters are k=43, l=10^6, err=0.12, etc.
        '''
        self.k = k
        self.l = l
        self.error_rate = err
        self.nu = ecc_len
        self.lbd = lbd
        self.xi = xi

        # ------------
        self.t = t
        self.ecc_msg_len = t + lbd * 2 + xi

        # ------------
        self.gamma = 1 - (err * 2)
        self.dec_check = ecc_len * (0.5 - self.gamma/3)

        # Irreducible polynomial for GF(2^128)
        self.irreducible_poly = galois.primitive_poly(2, 128)

        # SHA3-512 and corresponding parameters
        self.hash = hashlib.sha3_512()
        self.L = ceil((self.hash.digest_size * 8) / self.lbd) + 4
        self.hash = hashlib.sha3_512().name
        self.bigctxt = ''

        # File prefix
        self.file_prefix = file_prefix

        # Keep track of precise computation time (without disk IO)
        self.gen_timer = 0
        self.rep_timer = 0

        print("Done initializing")


    def bitarr(self, i):
        r = randbits(i)
        return [int(d) for d in bin(r)[2:].zfill(i)]
    
    def read_matrix(self, index):
        return np.load(f"LPN_Matrices/{index}.npy")

    # def LPN_enc_batch(self, A)
    # * precompute the noisy msg X self.l in GEN
    # * call the batch enc function with all the LPN matrices & the subsamples
    def LPN_batch_enc(self, keys, msgs):
        # Prep the messages to encode (put each on a new line)
        with open(f'src{self.file_prefix}.src', 'w') as f:
            f.writelines(msgs[0])
        
        # Encode message
        check_output(["./encode", "parity.pchk", "gen.gen", f"src{self.file_prefix}.src", f"e{self.file_prefix}.enc"])
        
        # Setup for adding errors to l encodings 
        code = []
        with open(f'e{self.file_prefix}.enc', 'r') as f:
            code = f.read().strip()
            code = [code+'\n'] * self.l
        
        with open(f'e{self.file_prefix}.enc', 'w') as f:
            f.writelines(code)

        # Adding errors
        # check_output(["./transmit", "e.enc", "r.rec", f"{randbits(32)}", "bsc", f"{self.error_rate}"])
        check_output(["bash", "ldpc_enc_batch.bash", f"-s {randbits(32)}", f"-e {self.error_rate}", "-f", self.file_prefix])

        # Reading noisy codes
        with open(f'r{self.file_prefix}.rec', 'r') as f:
            noisy_msg = f.readlines()

        # Transform the str output to a binary vector
        m = [np.array([int(b) for b in nm.strip()]) for nm in noisy_msg]
        
        # Multiply LPN matrices by the LPN keys (subsamples of iris code)
        t = time.time()
        ctxt = mx_par.gen_helper(self.l, keys, m)
        print(f"Computed {len(ctxt)} ciphertexts in {time.time() - t} seconds")

        # Compute l ciphetexts
        # ctxt = [m[i] ^ d[i] for i in range(self.l)]

        d = []

        return ctxt


    def LPN_dec_batch(self, As, keys, ctxts, process):
        '''
        Decypts a batch of ciphertexts in one go (in one call)
        '''

        # multiply LPN matrix by the key
        d = [np.matmul(self.read_matrix(As[i]), keys[i]) % 2 for i in range(len(ctxts))]

        # Add (mod 2) d and ctxt (assuming ctxt is a numpy array)
        tmp = ''
        for i in range(len(ctxts)):
            for j in ((d[i] ^ ctxts[i]) % 2):
                tmp += str(j)

        # print(f'Testing LPN_dec_batch. Process id: {process}\n temp: {len(tmp)}')
        
        # Encode temp into a bitstring
        input_file_name = f'r{self.file_prefix}-{process}.rec'
        with open(input_file_name, 'w') as f:
            f.write(tmp)
        
        # Call LDPC code (C) and decode temp
        decoded = check_output(["bash", "ldpc_dec_batch.bash", f"-c {input_file_name}", f"-e {self.error_rate}", "-p", str(process), "-f", self.file_prefix])
        decoded = decoded.decode("ASCII").split("\n")
        details, test = [i.split() for i in decoded[1:-3]], decoded[-3]
        # print(f'Testing LPN_dec_batch. Process id: {process}\n bash script output:',details)
        if (test.split()[3]) == '0':
            # print(f'Testing LPN_dec_batch. Process id: {process}. All invalid...')
            return np.array([]) # Invalid decryption
        
        print(f'Testing LPN_dec_batch. Process id: {process}. Found valid decoding!')
        valid_msg_index = -1
        for i in details:
            # Layout of i: [index_of_msg, num_of_iterations, valid/invalid (0/1), change%]; 
            if i[2] == '1':
                valid_msg_index = int(i[0])
                break
        

        output_file_name = f'e{self.file_prefix}-{process}.ext'
        with open(output_file_name, 'r') as f:
            out = f.readlines()

        decoded = out[valid_msg_index].strip()
        
        g_i = np.array([int(b) for b in decoded])
        # Check if hamming weigth of g_i is more than some value (depends on self.error_rate)
        # if sum(temp ^ g_i) > self.dec_check: # TODO!!!! FIGURE OUT HOW SUM WILL WORK WITH LIST OF ARRAYS
        #     return np.array([]) # This is our "error" output
        return g_i
    
    

    def mac(self, key, ciphertexts):
        t = time.time()
        h = hashlib.new(self.hash)

        if not self.bigctxt: self.bigctxt = str(ciphertexts)
        bigctxt = self.bigctxt

        # Generate a digest of entire ciphertexts
        bigctxt = bigctxt.encode()
        h.update(bigctxt)
        bigctxt = h.hexdigest()
        bigctxt = bin(int(bigctxt, base=16))[2:].zfill(512)

        # Encode digest into vector m
        m = [galois.Poly.Int(int(bigctxt[i:i+self.lbd], base=2)) for i in range(0, len(bigctxt), self.lbd)]

        # split key into x and y
        x = galois.Poly.Int(int(key[:self.lbd], base=2))
        y = galois.Poly.Int(int(key[self.lbd:], base=2))

        mx = mx_par.mx_serial(m, x, self.L-4, self.irreducible_poly) 

        T_rep = (pow(x, self.L, self.irreducible_poly) + (pow(x, 2, self.irreducible_poly) * mx) + (x * y)) % self.irreducible_poly
        
        print(f'Calculated MAC in {time.time() - t} seconds')

        return T_rep
        


    # w is the input string (iris encoding), n is the iris mask (used for sampling unmasked bits)
    def gen(self, w, n):
        # step 1: generate R and R_1
        R = bin(randbits(self.xi))[2:].zfill(self.xi)
        R_1 = bin(randbits(self.lbd * 2))[2:].zfill(self.lbd * 2)

        to_enc = (R + R_1).zfill(self.ecc_msg_len)

        # Generate l sets of unmasked positions 
        t = time.time()
        self.positions = [ random.SystemRandom().sample(np.flatnonzero(n).tolist(), k=self.k) for _ in range(self.l) ]
        print(f"Generating sample positions: {time.time() - t} sec ")
        samples = []
        # step 2: start a loop
        for i in range(self.l):
            # Get a sample of w at positions
            sample_i = np.array([w[pos] for pos in self.positions[i]])
            samples.append(sample_i)


        self.ctexts = self.LPN_batch_enc(samples, [to_enc for _ in range(self.l)])

        self.T = self.mac(R_1, self.ctexts) 

        # step 7: Output key R, self.ctexts, and self.T
        return R
        
    # w is W' in the paper, ciphertext and T can be found in self.ctexts and self.T respectively
    def rep_parallel(self, w, num_processes=1):
        # Pre-compute hash of ctxt TODO

        finished = multiprocessing.Array('b', False)
        a = np.array_split(range(self.l), 1000)
        b = np.array_split(range(1000), num_processes)
        finished = multiprocessing.Manager().list([None for _ in range(num_processes)])
        processes = []
        for x in range(num_processes):
            p = multiprocessing.Process(
                target=self.rep_process, args=(w, [a[i] for i in b[x]], finished, x)
            )
            processes.append(p)
            p.start()
        for p1 in processes:
            p1.join()
        if any(finished):
            print("Rep succeeded")
            return next(item for item in finished if item is not None)
        print("Rep failed")
        return None
    
    def rep_process(self, w_, arr_of_indices, finished, process_id):
        for indices in arr_of_indices:
            if any(finished):
                # print("One of the other threads returned")
                return 

            samples = []
            matrices = []
            ctxts = []
            for i in indices:
                sample_i = np.array([w_[pos] for pos in self.positions[i]])
                samples.append(sample_i)
                matrices.append(i)
                ctxts.append(self.ctexts[i])
            
            dec = self.LPN_dec_batch(matrices, samples, ctxts, process_id)

            if len(dec) > 0: print(dec[:15])

            # STEP iv
            if not (len(dec) == 0 or dec[:self.t].any()): # i.e., if dec is not None
                R = ''
                for c in dec[self.t:self.t + self.xi]:
                    R += str(c)
                
                R_1 = ''
                for c in dec[self.t + self.xi:]:
                    R_1 += str(c)

                T_rep = self.mac(R_1, self.ctexts)

                if T_rep == self.T:
                    print("Check passed")
                    finished[process_id] = R
                    return
        
        return



def img_opener(path, mask=False):
    ''' 
    INPUTS:
    `path`: String - Path to an image to be opened
    `mask`: Boolean - Optional flag to indicate whether the image is a mask
    OUTPUT: Flattened (1D) NumPy array of bits.
    Depending on `mask`, this is either just one or six array(s) of length 512 x 64
    '''
    image = Image.open(path)
    data = np.asarray(image) % 2 # map everything to 0 or 1 since it's B&W values (0 or 255)
    if mask:
        return data.flatten()
    else:
        return [data[i:i+64].flatten() for i in range(0, 384, 64)]


def main(first, toTest):
    mask1 = f"./NEWOutput/NormalizedMasks/{first}_mano.bmp"
    code1 = f"./NEWOutput/IrisCodes/{first}_code.bmp"

    m1 = img_opener(mask1, mask=True)
    c1 = [ m1 & c for c in img_opener(code1) ] # XOR all 6 codes (one per Gabor filter pair) with mask here


    print("Testing ", first)
    t1 = time.time()
    fe = FuzzyExtractor(l=1000000, file_prefix=first)
    t2 = time.time()
    print(f"Initialized (generated lpn arrays & GF(2^128)) in {t2 - t1} seconds")

    # HACK ACCORDING TO FULLERS PAPER (SECTION 4), TRANSFORM #5 HAS THE BEST RATE FOR IMAGES OF SAME IRIS
    a = fe.gen(c1[5], m1)
    t3 = time.time()
    print(f"Ran GEN in {t3 - t2} seconds")

    results = []

    for t in toTest:
        maskt = f"./NEWOutput/NormalizedMasks/{t}_mano.bmp"
        codet = f"./NEWOutput/IrisCodes/{t}_code.bmp"

        mt = img_opener(maskt, mask=True)
        ct = [ mt & c for c in img_opener(codet) ] # XOR all 6 codes (one per Gabor filter pair) with mask here

        t_ = time.time()
        b = fe.rep_parallel(ct[5], num_processes=multiprocessing.cpu_count())
        results.append(b)
        t1 = time.time()
        print(f"Ran REP parallel in {t1 - t_} seconds")



    print(a)
    print(results)


    print("no problems so far")


if __name__ == '__main__':
    sec = [
        '04841d617',
        '04841d616',
        '04841d613',
        '04841d610',
        '04841d305',
        '04841d312',
        '04841d591',
        '04841d586',
        '04841d602',
        '04841d606',
        '04841d603',
        '04841d307',
        '04841d320',
        '04841d308',
        '04841d590',
        '04841d313',
        '04841d593',
        '04841d594',
        '04841d597',
        '04841d314',
        '04841d321',
        '04841d601',
        '04841d589',
        '04841d604',
        '04841d319',
        '04841d612',
        '04841d306',
        '04841d311',
        '04841d618',
        '04841d588',
        '04841d596',
        '04841d595',
        '04841d309',
        '04841d310',
        '04841d605',
        '04841d607',
        '04841d315',
        '04841d600',
        '04841d592',
        '04841d611',
        '04841d587',
        '04841d599',
        '04841d318',
        '04841d615',
        '04841d598',
        '04841d317',
        '04841d609',
        '04841d316',
        '04841d614',
        '04841d608',
        '04841d304',
        ]
    main(first='04841d596', toTest=sec)

