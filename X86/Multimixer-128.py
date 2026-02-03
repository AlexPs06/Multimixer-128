# -*- coding: utf-8 -*-

# The following part is used to generate key and message

import random

def rand_key(p):
   
    # Variable to store the string
    key1 = ""
 
    # Loop to find the string of desired length
    for i in range(p):
         
        # randint function to generate 0, 1 randomly and converting the result into str
        tem = str(random.randint(0,1))
 
        # Concatenation the random 0, 1 to the final result
        key1 += tem
         
    return(key1)
 
#This is the block function for Multimixer-128

def Int_multimix(M,K):
    # Define the number of blocks and then Split the message and the key into blocks of size 256
    # Store the split msg blocks and key blocks in the array Msg_blks, Key_blks respectively
    Block_num = len(M)//256
    Msg_blks = []
    Key_blks = []
    for i in range(Block_num):
        Msg_blk = M[256*i: 256*(i+1)]
        Key_blk = K[256*i: 256*(i+1)]
        Msg_blks.append(Msg_blk)
        Key_blks.append(Key_blk)     
    
    # For each block, divide the messages and keys into 32 bit chunks to obtain the A,B,P and Q
    Blk_reslts = []
    for i in range(len(Msg_blks)):
        # X_i's are the first parts of the messages, Y_i's are the second part, similar for H_i and K_i. 
        
        X_0 = int(Msg_blks[i][0:32],2)
        X_1 = int(Msg_blks[i][32:64],2)
        X_2 = int(Msg_blks[i][64:96],2)
        X_3 = int(Msg_blks[i][96:128],2)

        Y_0 = int(Msg_blks[i][128:160],2)
        Y_1 = int(Msg_blks[i][160:192],2)
        Y_2 = int(Msg_blks[i][192:224],2)
        Y_3 = int(Msg_blks[i][224:256],2)

        H_0 = int(Key_blks[i][0:32],2)
        H_1 = int(Key_blks[i][32:64],2)
        H_2 = int(Key_blks[i][64:96],2)
        H_3 = int(Key_blks[i][96:128],2)
        
        K_0 = int(Key_blks[i][128:160],2)
        K_1 = int(Key_blks[i][160:192],2)
        K_2 = int(Key_blks[i][192:224],2)
        K_3 = int(Key_blks[i][224:256],2)
        #A_i and B_i are the inputs to the block function (after addition with the corresponding key bits)

        # print(hex(X_0))
        # print(hex(X_1))
        # print(hex(X_2))
        # print(hex(X_3))

        # print(Y_0)
        # print(Y_1)
        # print(Y_2)
        # print(Y_3)
        
        # print(H_0)
        # print(H_1)
        # print(H_2)
        # print(H_3)

        # print(hex(H_0))
        # print(hex(H_1))
        # print(hex(H_2))
        # print(hex(H_3))

        # print(K_0)
        # print(K_1)
        # print(K_2)
        # print(K_3)

        A_0 =(X_0+H_0)%(2**32) #X
        A_1 =(X_1+H_1)%(2**32)
        A_2 =(X_2+H_2)%(2**32)
        A_3 =(X_3+H_3)%(2**32)
        
        

        B_0 =(Y_0+K_0)%(2**32) #Y
        B_1 =(Y_1+K_1)%(2**32)
        B_2 =(Y_2+K_2)%(2**32)
        B_3 =(Y_3+K_3)%(2**32)
        # Compute P_i and Q_i from A_i, B_i 
        P_0 = (A_0+A_1+A_2)%(2**32)
        P_1 = (A_1+A_2+A_3)%(2**32)
        P_2 = (A_2+A_3+A_0)%(2**32)
        P_3 = (A_3+A_0+A_1)%(2**32)

        

        Q_0 = (B_1+B_2+B_3)%(2**32)
        Q_1 = (B_2+B_3+B_0)%(2**32)
        Q_2 = (B_3+B_0+B_1)%(2**32)
        Q_3 = (B_0+B_1+B_2)%(2**32)


        # print(hex(A_0), hex(A_1), hex(A_2), hex(A_3))
        # print(hex(B_0), hex(B_1), hex(B_2), hex(B_3))

        # print(hex(Y_0), hex(Y_1), hex(Y_2), hex(Y_3))
        # print(hex(K_0), hex(K_1), hex(K_2), hex(K_3))
        # print(hex(H_0), hex(H_1), hex(H_2), hex(H_3))
        # print(hex(P_0), hex(P_1), hex(P_2), hex(P_3))
        # print(hex(Q_0), hex(Q_1), hex(Q_2), hex(Q_3))

        # the block function computes the 8 multiplications
        # Results for each block are stored as 8-tuple in the array Blk_reslts
        Blk_res = (A_0*B_0, A_1*B_1, A_2*B_2, A_3*B_3, P_0*Q_0, P_1*Q_1, P_2*Q_2, P_3*Q_3)
        
        # print(hex(A_0))
        # print(hex(B_0))
        # print(hex(A_0*B_0))

        # a = 0x36fe5b2d
        # b = 0xb30d1b3c
        # print(hex(a*b))
        # print(hex(A_1))
        # print(hex(B_1))
        # print(hex(A_1*B_1))

        # print(hex(A_2))
        # print(hex(B_2))
        # print(hex(A_2*B_2))
        # print(hex(A_0*B_0),hex(A_1*B_1),hex(A_2*B_2),hex(A_3*B_3),hex(P_0*Q_0),hex(P_1*Q_1),hex(P_2*Q_2),hex(P_3*Q_3))

        # exit(1);

        Blk_reslts.append(Blk_res)
    # res finally adds the results of each block co-ordinatewise and then converts it back to a string
    res =""
    for j in range(8):
        temp = '{:064b}'.format(sum(i[j] for i in Blk_reslts)%(2**(64)))
        res += temp
    
#print(Blk_reslts)

    return(res)


#256 = block size of Multimixer-128


# l = int(input("Enter message Length:   "))            

# Msg = rand_key(l*256)
# Key = rand_key(l*256)

# print(Msg)
# print(Key)
Msg = "0000101000101000101101101101000001111010111110000001101000010001011110011010000100100010101111011011011110000111010101010011100100110011010010000000010100101000000111010110001111100011010000010011101000100011101001011110011001110000101111010001100010011100"
Key = "1100101111001011010111010011001010111100000001100100000100011100100010100111011100101001111001011100110100011101100000101001011100100011000110011010110010011011100101011010100100110111111110110111001101001000001011101001110100110011101101101001001011001111"
print(Int_multimix(Msg, Key))