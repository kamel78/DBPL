# DBPL
Delphi’s Bilinear Pairings Library (DBPL)
DBPL is a Delphi components library that implements bilinear pairings functionalities on several pairing friendly curves. 

The developed library implements bilinear pairing functions on different curves using different possible sets of parameters. The first level of the library implements Multi-precision arithmetic over Fp using assembly instructions in both 32bit and 64bit modes to gain high runtime performances, while best optimal algorithms for squaring, division, modular inversion quadratic residues has been chosen.  In addition, arithmetic over extension fields is implemented over several extensions (Fp3,Fp4…,Fp36) is implemented optimally using best known algorithms. In a second level, arithmetic over elliptic curves is implemented with respect to the three coordinates systems   : affine, projective and Jacobean. On the top of the described routines, an exhaustive set of pairings friendly curves is implemented using several optimal parameters, including:
-	Super singular curves;
-	MNT curves
-	BN curves: including parameters presented by Beuchat, Aranha, Scott and the ISO normalized curves.
-	BLS curves: for the families    BLS12, BLS24 and BLS27.
-	KSS curves: including KSS16, KSS18 and KSS36.
For each curve, different possible pairing algorithms are implemented (Eta, Tate, Ate and Opt-Ate), with respect to several coordinates systems, and several predefined parameters. The final exponentiation step is optimized using Karabiner’s algorithm. The whole set of functions is implemented using oriented-object model as a set of Delphi’s non-visual components, in order to provide simple and assisted use of the library. This library provides common additionally a visual components inherited from TreeView that presents the parameters of each pairing friendly curve when linked to the corresponding non-visual component. Already supports two main platforms: 64 and 32 bit.
Please refer to DPBL-Demonstration for more examples (a little and simple demo is also provided in LittleDemo).
The library is compatible with Delphi Xe2-Xe8 and Delphi DX-10.1-10.2

Components :
![alt text](https://github.com/kamel78/DBPL/blob/master/Component.png)

Demo :
![alt text](https://github.com/kamel78/DBPL/blob/master/Appdemo.png)

Features
With DPBL, you can implement these features easily:
-	Generating curve’s point on Fields and Extensions-Fields by direct hashing from strings
-	Simple choice and initialization of pairings parameters (either from the object’s inspector or by code)
-	Pairings with different algorithms and according to the three possible coordinates systems (Affine, Projective and Jacobians)
-	Handle both 32bit and 64bit architectures (Native Pascal can also be automatically used for Android/IOS implementations)   
-	Extremely fast pairing….
-	Display of results (integer’s values, points ‘coordinates, extension fields values..) in decimal/hexadecimal forms as strings.
-	Support of customized parameterization  an initialization
-	
Usage
Can be used to implements different variants of Identity based cryptosystems, pairing based Signature schemes, IBS, ABS…… 
•	
Contribute
Welcome contribution! 

