修論　サイドチャネル攻撃に安全かつ高速な耐量子暗号CSIDH

Stronger and Faster Side-Channel Protections for CSIDH (https://github.com/JJChiDguez/csidh)

CSI-Fish Efficient isogeny based signatures through class group computations (https://github.com/KULeuven-COSIC/CSI-FiSh)

のソースを組み合わせて作成．

action_costについて，babaiのアルゴリズムなどはカウントしてない．

前半の計算において計算の仕方を変えることによって，後半の計算の同種写像計算の回数が変わる．VECTOR1とVECTOR2はそのやり方の違い．VECTOR2の方が若干早い．

# ------------------------------------------------------------------------
# C code
To compile the files you can do the following:

First, you can use any version of gcc compiler (just set the variable CC as 
an input of the Makefile [variable CC is optional, gcc is set by default]).

# Testing a CSIDH protocol (key exchange protocol)
[Compilation]

	(Using dummy operations and one torsion point)
		make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1
	(Using dummy operations and two torsion points)
		make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2
	(Dummy-free approach and using two torsion points)
		make csidh BITLENGTH_OF_P=512 TYPE=DUMMYFREE
	(提案手法１)
		make csidh BITLENGTH_OF_P=512 TYPE=PROPOSAL_1
	(提案手法２)
		make csidh BITLENGTH_OF_P=512 TYPE=PROPOSAL_2
		

[Execution]

		./bin/csidh


# Running-time: number of field operations
[Compilation]

	(Using dummy operations and one torsion point)
		make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1
	(Using dummy operations and two torsion points)
		make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2
	(Dummy-free approach and using two torsion points)
		make action_cost BITLENGTH_OF_P=512 TYPE=DUMMYFREE
	(提案手法１)
		make action_cost BITLENGTH_OF_P=512 TYPE=PROPOSAL_1
	(提案手法２)
		make action_cost BITLENGTH_OF_P=512 TYPE=PROPOSAL_2

[Execution]

		./bin/action_cost

# Running-time: number of clock cycles
[Compilation]

	(Using dummy operations and one torsion point)
		make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1
	(Using dummy operations and two torsion points)
		make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2
	(Dummy-free approach and using two torsion points)
		make action_timing BITLENGTH_OF_P=512 TYPE=DUMMYFREE
	(提案手法１)
		make action_timing BITLENGTH_OF_P=512 TYPE=PROPOSAL_1
	(提案手法２)
		make action_timing BITLENGTH_OF_P=512 TYPE=PROPOSAL_2


[Execution]

		./bin/action_timing

# Clean

	make clean
