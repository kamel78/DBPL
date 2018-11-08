unit Fp12Arithmetic;


interface

uses System.SysUtils,VCL.dialogs, LargeIntegers, Fp2Arithmetic,Fp6Arithmetic, GeneralTypes;

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp12. Fp12 is the Twelvtic extention of Fp (Tower extention of order 2 for FP6)
   with respect to the irreductible polynômial W^2-Gamma=0. Elements are in the Form a+b*W (a and b are from Fp6).

********************************************************************************************}

Type


 PFp12Int=^Fp12Int;

  Fp12Int=record
  a,b:Fp6Int;
  Tower:PtrTowerParams12;
  public
    class operator Add(const Left, Right: Fp12Int): Fp12Int;
    class operator Subtract(const Left, Right: Fp12Int): Fp12Int;
    class operator Multiply(const Left, Right: Fp12Int): Fp12Int;
    class operator Equal(const Left, Right: Fp12Int): Boolean;
    class operator NotEqual(const Left, Right: Fp12Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse:Fp12Int;
    function Sqr:Fp12Int;
    function Conjugate:Fp12Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function Pow(e:LInt;powerMode:FpPoweringMode=pmNormal):Fp12Int;
    function PowToP:Fp12Int;                      /// optimized power computation using Frobenius Map
    function PowToP2:Fp12Int;
    function PowToP3:Fp12Int;
    procedure SetTowerParams(param:PtrTowerParams12);
    procedure SetToOne;
    procedure SetToZero;
    procedure SetToRandom;
    procedure SetFromStrings(a1,b1,c1,a2,b2,c2:String);
  end;

    procedure _FP6_From_Strings(const valueA1,valueB1,valueC1,valueA2,valueB2,valueC2:String;var Result:Fp12Int);
    function _FP12_To_DecimalString(const value:Fp12Int):string;
    procedure _Inv_FP12(const value:Fp12Int;var Result:Fp12Int);
    procedure _Pow_FP12_P3(const value:Fp12Int;var Result:Fp12Int);
    procedure _Pow_FP12_P2(const value:Fp12Int;var Result:Fp12Int);
    procedure _Pow_FP12_P(const value:Fp12Int;var Result:Fp12Int);
    procedure _Pow_FP12(const value:Fp12Int;Exponent:LInt;var Result:Fp12Int;Poweringmode :FpPoweringMode=pmNormal);
    procedure _Conjugate_FP12(const Value:Fp12Int;var result:Fp12Int);
    procedure _Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
    procedure _Be_Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
    procedure _Mul_FP12(const Left,Right:Fp12Int; var Result:Fp12Int);
    procedure _Sub_FP12(const Left,Right:Fp12Int; var Result:Fp12Int);
    procedure _Add_FP12(const Left,Right:Fp12Int; var Result:Fp12Int);
    function _Equals_FP12(const Left,Right:Fp12Int):boolean;
    procedure _Compressed_Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
    procedure _DeCompressed_Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
    procedure _Sparse_Mul_FP12(const Left,Right:Fp12Int; var Result:Fp12Int;TwistMode:TTwistModel);
    procedure _Mul_FP12_By_W(const value: Fp12Int;Gamma:Fp6Int; var Result: Fp12Int);
    procedure _Neg_FP12(const Value: Fp12Int; var Result: Fp12Int);

implementation

{*******************************************************************************}
            ///      Procedures for FP12 Arithmetic
{*******************************************************************************}

        {**********   Add two FP12 integers *****************}
procedure _Add_FP12(const Left,Right:Fp12Int; var Result:Fp12Int);
begin
Result.Tower:=Left.Tower;
_Add_FP6(Left.a,Right.a,Result.a);
_Add_FP6(Left.b,Right.b,Result.b);
end;

        {**********   Sub two FP12 integers *****************}
procedure _Sub_FP12(const Left,Right:Fp12Int; var Result:Fp12Int);
begin
Result.Tower:=Left.Tower;
_Sub_FP6(Left.a,Right.a,Result.a);
_Sub_FP6(Left.b,Right.b,Result.b);
end;

        {********** Multiply two FP12 integers *****************}
procedure _Sparse_Mul_FP12(const Left,Right:Fp12Int; var Result:Fp12Int; TwistMode:TTwistModel);
var t:array[0..4] of Fp6Int;
    Sigmaisiplus1:boolean;
begin
Sigmaisiplus1:=(Left.Tower.Sigma.a.Data.i16[-2]=1)and(LEft.Tower.Sigma.b.Data.i16[-2]=1)and(Left.Tower.Sigma.a.Data.i16[0]=1)and(Left.Tower.Sigma.b.Data.i16[0]=1);
Result.Tower:=Left.Tower;
if TwistMode=twDType then begin
                          _Sparse1_Mul_FP6(Left.a,Right.a,t[0]);
                          _Sparse2_Mul_FP6(Left.b,Right.b,t[1]);
                          end
else begin
     _Sparse2_Mul_FP6(Left.a,Right.a,t[0]);
     _Sparse3_Mul_FP6(Left.b,Right.b,t[1]);
     end;
_Add_FP6(Left.a,Left.b,t[2]);
if TwistMode=twDType then _Sparse_Add_FP6(Right.a,Right.b,t[3])
else _Add_FP6(Right.b,Right.a,t[3]);
if Sigmaisiplus1 then begin
                      /// Works Only if Sigma is 1+i
                      t[4].Tower:=t[1].Tower;
                      t[4].b:=t[1].a;
                      t[4].c:=t[1].b;
                      _Sub_Lint(t[1].c.a,t[1].c.b,t[4].a.a);
                      _Add_Lint(t[1].c.a,t[1].c.b,t[4].a.b);
                      ///
                      end
else _Mul_FP6_By_V2(t[1],t[4]);
_Add_FP6(t[0],t[4],Result.a);
_Sparse2_Mul_FP6(t[2],t[3],Result.b);
_Nomod_Sub_FP6(Result.b,t[0],Result.b);
_Sub_FP6(Result.b,t[1],Result.b);
end;

        {********** Multiply two FP12 integers *****************}
procedure _Mul_FP12(const Left,Right:Fp12Int; var Result:Fp12Int);
var t:array[0..4] of Fp6Int;
begin
Result.Tower:=Left.Tower;
_Mul_FP6(Left.a,Right.a,t[0]);
_Mul_FP6(Left.b,Right.b,t[1]);
_Add_FP6(Left.a,Left.b,t[2]);
_Add_FP6(Right.a,Right.b,t[3]);
_Mul_FP6_By_V2(t[1],t[4]);
_Add_FP6(t[0],t[4],Result.a);
_Mul_FP6(t[2],t[3],Result.b);
_Sub_FP6(Result.b,t[0],Result.b);
_Sub_FP6(Result.b,t[1],Result.b);
end;

        {********** Multiply FP12 with W *****************}
procedure _Mul_FP12_By_W(const value: Fp12Int;Gamma:Fp6Int; var Result: Fp12Int);  // Value should be different than Result
begin
Result.Tower:=Value.Tower;
Result.b:=Value.a;
_Mul_FP6(Value.b,Gamma,Result.a);
end;

        {********** Get Square of an FP12 *****************}
        // Classical approach: for Square and multiply algorithm
procedure _Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
var  t:array[0..1] of Fp6Int;
      Sigmaisiplus1:boolean;
begin
      {       Using Karatsuba Squering    }
Result.SetTowerParams(Value.Tower);
{
_Sqr_FP6(Value.a,t[0]);
_Sqr_FP6(Value.b,t[1]);
//_Mul_FP6(value.a,Value.b,t[2]);
_Mul_FP6(value.a,Value.b,Result.b);
_Mul_FP6_Gamma(t[1],Result.a);
_Add_FP6(Result.a,t[0],Result.a);
//_Add_FP6(t[2],t[2],Result.b);
_Shl_LInt(Result.b.a.a,1);
_Shl_LInt(Result.b.a.b,1);
_Shl_LInt(Result.b.b.a,1);
_Shl_LInt(Result.b.b.b,1);
_Shl_LInt(Result.b.c.a,1);
_Shl_LInt(Result.b.c.b,1);
_Mod_LInt(Result.b.a.a, Result.Tower.FieldParam.p,Result.b.a.a);
_Mod_LInt(Result.b.a.b,Result.Tower.FieldParam.p,Result.b.a.b);
_Mod_LInt(Result.b.b.a,Result.Tower.FieldParam.p,Result.b.b.a);
_Mod_LInt(Result.b.b.b,Result.Tower.FieldParam.p,Result.b.b.b);
_Mod_LInt(Result.b.c.a,Result.Tower.FieldParam.p,Result.b.c.a);
_Mod_LInt(Result.b.c.b,Result.Tower.FieldParam.p,Result.b.c.b); }

// This is a Faster way to sequare FP12 that optimise if  Sigma=1+i (Generally the case)!
Sigmaisiplus1:=(Value.Tower.Sigma.a.Data.i16[-2]=1)and(Value.Tower.Sigma.b.Data.i16[-2]=1)and(Value.Tower.Sigma.a.Data.i16[0]=1)and(Value.Tower.Sigma.b.Data.i16[0]=1)and(Value.Tower.FieldParam.Beta=-1);
_Add_FP6(value.a,Value.b,t[0]);
t[1].SetTowerParams(value.Tower);
if Sigmaisiplus1 then begin
                      _Sub_Lint(Value.b.c.a,Value.b.c.b,t[1].a.a);
                      _Add_Lint(Value.b.c.a,Value.b.c.b,t[1].a.b);
                      end
else _Mul_FP6_By_V2(Value.b,t[1]);
_Add_FP2(t[1].a,Value.a.a,t[1].a);
_Add_FP2(value.b.a,value.a.b,t[1].b);
_Add_FP2(Value.b.b,Value.a.c,t[1].c);
_Mul_FP6(value.a,value.b,result.b);
_Mul_FP6(t[0],t[1],Result.a);
if Sigmaisiplus1 then begin
                      _Sub_Lint(Result.b.c.a,Result.b.c.b,t[1].a.a);
                      _Add_Lint(Result.b.c.a,Result.b.c.b,t[1].a.b);
                      end
else _Mul_FP6_By_V2(Result.b,t[1]);
_Add_FP2(t[1].a,Result.b.a,t[1].a);
_Add_FP2(Result.b.a,Result.b.b,t[1].b);
_Add_FP2(Result.b.b,Result.b.c,t[1].c);
_Sub_FP6(Result.a,t[1],Result.a);
_Add_FP6(Result.b,Result.b,Result.b);
end;

        {********** Get Square of an FP4 element *****************}
        // Alternative Defined for Beuchat's Algorithm ovec Cyclotomic Group
procedure _Sqr_Fp4(const a0,a1:Fp2Int;tower:PtrTowerParams12; var c0,c1:Fp2Int);
var t:array[0..1] of Fp2Int;
begin
_Sqr_FP2(a0,t[0]);
_Sqr_FP2(a1,t[1]);
_Mul_FP2(t[1],tower^.Sigma,c0);
_Add_FP2(c0,t[0],c0);
_Add_FP2(a0,a1,c1);
_Sqr_FP2(c1,c1);
_Sub_FP2(c1,t[0],c1);
_Sub_FP2(c1,t[1],c1);
end;

     {********** Get Compressed Square of an FP12 (Karabina Approach)*****************}
     // Value should verify the condition (Value)^(p^6+1)=1 (Cyclotomic)
procedure _Compressed_Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
var  t:array[0..6] of Fp2Int;
begin
_Sqr_FP2(value.a.b,t[0]);
_Sqr_FP2(value.b.c,t[1]);
_Add_FP2(value.a.b,value.b.c,t[5]);
_Sqr_FP2(t[5],t[2]);
_Nomod_Add_FP2(t[0],t[1],t[3]);
_Nomod_Sub_FP2(t[2],t[3],t[3]);
_Mod_FP2_FP(t[3],value.Tower.FieldParam.p,t[5]);
_Add_FP2(value.b.a,value.a.c,t[6]);
_Sqr_FP2(t[6],t[3]);
_Sqr_Fp2(value.b.a,t[2]);
_Mul_FP2(t[5],value.Tower.Sigma,t[6]);
_Nomod_Add_FP2(t[6],value.b.a,t[5]);
_Nomod_Add_FP2(t[5],t[5],t[5]);
_Add_FP2(t[5],t[6],Result.b.a);
_Mul_FP2(t[1],value.Tower.Sigma,t[4]);
_Nomod_Add_FP2(t[0],t[4],t[5]);
_Nomod_Sub_FP2(t[5],value.a.c,t[6]);
_Add_FP2(t[6],t[6],t[6]);
_Sqr_FP2(value.a.c,t[1]);
_Add_FP2(t[5],t[6],result.a.c);
_Mul_FP2(t[1],value.Tower.Sigma,t[4]);
_Add_FP2(t[4],t[2],t[5]);
_Nomod_Sub_FP2(t[5],value.a.b,t[6]);
_Add_FP2(t[6],t[6],t[6]);
_Add_FP2(t[5],t[6],result.a.b);
_Nomod_Add_FP2(t[1],t[2],t[0]);
_Nomod_Sub_FP2(t[3],t[0],t[5]);
_Nomod_Add_FP2(t[5],value.b.c,t[6]);
_Nomod_Add_FP2(t[6],t[6],t[6]);
_Add_FP2(t[6],t[5],Result.b.c);
Result.b.b.a:=0;
Result.b.b.b:=0;
Result.a.a.a:=0;
Result.a.a.b:=0;
end;

    {********** Get Decompressed Square of an FP12 (Karabina Approach)*****************}
procedure _DeCompressed_Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
var t:array[0..6] of Fp2Int;
begin
if _Is_FP2_Null(value.a.b) then begin
                                _Mul_FP2(value.a.b,value.b.c,t[1]);
                                _Add_FP2(t[1],t[1],t[1]);
                                _Inv_FP2(value.a.c,t[2]);
                                _Mul_FP2(t[1],t[2],result.b.b);
                                _Sqr_FP2(t[1],t[2]);
                                _Add_FP2(t[2],t[2],t[2]);
                                _Mul_FP2(value.a.c,value.a.b,t[3]);
                                _Add_FP2(t[3],t[3],t[4]);
                                _Add_FP2(t[3],t[4],t[4]);
                                _Sub_FP2(t[3],t[4],t[0]);
                                _Mul_FP2(t[0],value.Tower.Sigma,Result.a.a);
                                _Inc_LInt(Result.a.a.a,1);
                                end
else begin
     _Sqr_FP2(value.b.c,t[5]);
     _Mul_FP2(t[5],value.Tower.Sigma,t[5]);
     _Sqr_FP2(value.a.b,t[4]);
     _Add_FP2(t[4],t[4],t[6]);
     _Add_FP2(t[4],t[6],t[4]);
     _Add_FP2(value.a.c,value.a.c,t[3]);
     _Add_FP2(value.b.a,value.b.a,t[2]);
     _Add_FP2(t[2],t[2],t[2]);
     _Inv_FP2(t[2],t[2]);
     _Add_FP2(t[4],t[5],t[1]);
     _Sub_FP2(t[1],t[3],t[1]);
     _Mul_FP2(t[1],t[2],Result.b.b);
     _Sqr_FP2(Result.b.b,t[1]);
     _Add_FP2(t[1],t[1],t[1]);
     _Mul_FP2(value.b.a,value.b.c,t[2]);
     _Mul_FP2(value.a.c,value.a.b,t[3]);
     _Add_FP2(t[3],t[3],t[4]);
     _Add_FP2(t[3],t[4],t[3]);
     _Add_FP2(t[1],t[2],t[0]);
     _Sub_FP2(t[0],t[3],t[0]);
     _Mul_FP2(t[0],value.Tower.Sigma,Result.a.a);
     _Inc_LInt(Result.a.a.a,1);
     end;
Result.a.b:=Value.a.b;
Result.a.c:=Value.a.c;
Result.b.a:=Value.b.a;
Result.a.b:=Value.a.b;
Result.b.c:=Value.b.c;
end;

        {********** Get Square of an FP12 (Beuchat Approach)*****************}
          // Value should verify the condition (Value)^(p^6+1)=1
procedure _Be_Sqr_FP12(const value: Fp12Int; var Result: Fp12Int);
var  t00,t11,t10,t01,t02,t12,aux,tmp: Fp2Int;
begin
      {       Using Beuchat FP4 technique    }
Result.Tower:=Value.Tower;
_Sqr_Fp4(value.a.a,value.b.b,value.Tower,t00,t11);
_Sqr_Fp4(value.b.a,value.a.c,value.Tower,t01,t12);
_Sqr_Fp4(value.a.b,value.b.c,value.Tower,t02,aux);
_Mul_FP2(aux,value.Tower^.Sigma,t10);
_Add_FP2(t00,t00,tmp);
_Add_FP2(t00,tmp,tmp);
_Sub_FP2(tmp,value.a.a,tmp);
_Sub_FP2(tmp,value.a.a,Result.a.a);
_Add_FP2(t01,t01,tmp);
_Add_FP2(t01,tmp,tmp);
_Sub_FP2(tmp,value.a.b,tmp);
_Sub_FP2(tmp,value.a.b,Result.a.b);
_Add_FP2(t02,t02,tmp);
_Add_FP2(t02,tmp,tmp);
_Sub_FP2(tmp,value.a.c,tmp);
_Sub_FP2(tmp,value.a.c,Result.a.c);
_Add_FP2(t10,t10,tmp);
_Add_FP2(t10,tmp,tmp);
_Add_FP2(tmp,value.b.a,tmp);
_Add_FP2(tmp,value.b.a,Result.b.a);
_Add_FP2(t11,t11,tmp);
_Add_FP2(t11,tmp,tmp);
_Add_FP2(tmp,value.b.b,tmp);
_Add_FP2(tmp,value.b.b,Result.b.b);
_Add_FP2(t12,t12,tmp);
_Add_FP2(t12,tmp,tmp);
_Add_FP2(tmp,value.b.c,tmp);
_Add_FP2(tmp,value.b.c,Result.b.c);
end;

        {**********   Conjugate an FP12 integer *************}
procedure _Conjugate_FP12(const Value:Fp12Int;var result:Fp12Int);
begin
Result.Tower:=Value.Tower;
Result.a:=Value.a;
_Neg_FP6(Value.b,Result.b);
end;

        {********** Raise an FP12 to a LInt power ************}
procedure _Pow_FP12(const value:Fp12Int;Exponent:LInt;var Result:Fp12Int;Poweringmode :FpPoweringMode=pmNormal);
var i:word;
    zi,tmp:Fp12Int;
    nafexpo,t1,t2:LIntArrayForm;
    SizeNaf,Sizebin:Integer;
begin
Result.Tower:=Value.Tower;
case Poweringmode of
pmNormal:begin
              //  Raise an FP12 to a LInt power (Classical Square and Multiply approach)************
         Result:=Value;
         for i:=Exponent.BitLength-2 downto 0 do begin
                                                 _Sqr_FP12(Result,Result);
                                                 if _Is_BitSet_At(Exponent,i) then _Mul_FP12(Result,value,Result);
                                                 end;
         end;
pmBeuchat:begin
                  //  Raise an FP12 to a LInt power (Beucaht FP4 Approach)************
                  //  Value should verify the condition (Value)^(p^6+1)=1
         Result:=Value;
         for i:=Exponent.BitLength-2 downto 0 do begin
                                                 _Be_Sqr_FP12(Result,Result);
                                                 if _Is_BitSet_At(Exponent,i) then _Mul_FP12(Result,value,Result);
                                                 end;
          end;
pmKarbina:begin
                  //  Raise an FP12 to a LInt power (Karabina Approach)************
                  //  Value should verify the condition (Value)^(p^6+1)=1
          zi.SetTowerParams(Value.Tower);
          tmp:=value;
          t1:=LIntToNAF(Exponent);
          t2:=LIntToIntArray(Exponent);
          SizeNaf:=0;
          SizeBin:=0;
          for i:=0 to length(t1)-1 do if t1[i]<>0 then inc(SizeNaf);
          for i:=0 to length(t2)-1 do if t2[i]<>0 then inc(SizeBin);
          if SizeNaf>SizeBin then NafExpo:=t2
          else NafExpo:=t1;
          if NafExpo[0]<>0  then begin
                                     Result:=Value;
                                     if NafExpo[0]=-1 then _Conjugate_FP12(Result,Result);
                                     end
          else begin
               Result.SetToOne;
               Result.SetTowerParams(Value.Tower);
               end;
          for i:=1 to Length(NafExpo)-1 do begin
                                           _Compressed_Sqr_FP12(tmp,tmp);
                                           if nafexpo[i]=1 then begin
                                                                _DeCompressed_Sqr_FP12(tmp,zi);
                                                                _Mul_FP12(zi,Result,Result);
                                                                end
                                           else if nafexpo[i]=-1 then begin
                                                                      _DeCompressed_Sqr_FP12(tmp,zi);
                                                                      _Conjugate_FP12(zi,zi);
                                                                      _Mul_FP12(zi,Result,Result);
                                                                      end
                                           end;
          end;
end;
if _IsNeg(Exponent) then _Conjugate_FP12(Result,Result);
end;

        {********** Raise an FP12 to " p " power ************}
        /// use Precomputed Frobenius Constants
procedure _Pow_FP12_P(const value:Fp12Int;var Result:Fp12Int);
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.a.a:=Value.a.a.Conjugate;
_Mul_FP2(Value.a.b.Conjugate,Value.Tower^.FrobeniusP_Const[1],Result.a.b);
_Mul_FP2(Value.a.c.Conjugate,Value.Tower^.FrobeniusP_Const[3],Result.a.c);
_Mul_FP2(Value.b.a.Conjugate,Value.Tower^.FrobeniusP_Const[0],Result.b.a);
_Mul_FP2(Value.b.b.Conjugate,Value.Tower^.FrobeniusP_Const[2],Result.b.b);
_Mul_FP2(Value.b.c.Conjugate,Value.Tower^.FrobeniusP_Const[4],Result.b.c);
end;

        {********** Raise an FP12 to " p^2 " power ************}
        /// use Precomputed Frobenius Constants
procedure _Pow_FP12_P2(const value:Fp12Int;var Result:Fp12Int);
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.a.a:=Value.a.a;
_Mul_FP2(Value.a.b,Value.Tower^.FrobeniusP2_Const[1],Result.a.b);
_Mul_FP2(Value.a.c,Value.Tower^.FrobeniusP2_Const[3],Result.a.c);
_Mul_FP2(Value.b.a,Value.Tower^.FrobeniusP2_Const[0],Result.b.a);
_Mul_FP2(Value.b.b,Value.Tower^.FrobeniusP2_Const[2],Result.b.b);
_Mul_FP2(Value.b.c,Value.Tower^.FrobeniusP2_Const[4],Result.b.c);
end;

        {********** Raise an FP12 to " p^3 " power ************}
        /// use Precomputed Frobenius Constants
procedure _Pow_FP12_P3(const value:Fp12Int;var Result:Fp12Int);
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.a.a:=Value.a.a.Conjugate;
_Mul_FP2(Value.a.b.Conjugate,Value.Tower^.FrobeniusP3_Const[1],Result.a.b);
_Mul_FP2(Value.a.c.Conjugate,Value.Tower^.FrobeniusP3_Const[3],Result.a.c);
_Mul_FP2(Value.b.a.Conjugate,Value.Tower^.FrobeniusP3_Const[0],Result.b.a);
_Mul_FP2(Value.b.b.Conjugate,Value.Tower^.FrobeniusP3_Const[2],Result.b.b);
_Mul_FP2(Value.b.c.Conjugate,Value.Tower^.FrobeniusP3_Const[4],Result.b.c);
end;

        {********** Get inverse of an FP12 *****************}
procedure _Inv_FP12(const value:Fp12Int;var Result:Fp12Int);
var t:array[0..2]of Fp6Int;
begin
Result.Tower:=Value.Tower;
_Sqr_FP6(value.a,t[0]);
_Sqr_FP6(value.b,t[1]);
_Mul_FP6_By_V2(t[1],t[2]);
_Sub_FP6(t[0],t[2],t[0]);
_Inv_FP6(t[0],t[1]);
_Mul_FP6(Value.a,t[1],Result.a);
_Mul_FP6(Value.b,t[1],Result.b);
_Neg_FP6(Result.b,Result.b);
end;

        {********** Get negative of an FP12 *****************}
procedure _Neg_FP12(const Value: Fp12Int; var Result: Fp12Int);
begin
Result.Tower:=Value.Tower;
_Neg_FP6(Value.a,Result.a);
_Neg_FP6(Value.b,Result.b);
end;

        {**********   Compare two FP12 integers *************}
function _Equals_FP12(const Left,Right:Fp12Int):boolean;
begin
result:=_Equals_FP6(Left.a,Right.a) and _Equals_FP6(Left.b,Right.b);
end;

        {****** Convert an FP12 to a Decimal String **********}
function _FP12_To_DecimalString(const value:Fp12Int):string;
begin
if Value.a.IsZero then begin
                       if value.b.IsZero then result:='0'
                       else result:='('+value.b.ToDecimalString+')* w';
                       end
else begin
     if value.b.IsZero then result:=value.a.ToDecimalString
     else Result:=Value.a.ToDecimalString+'+('+value.b.ToDecimalString+')* w';
     end;


{Result:=Value.a.a.a.ToDecimalString+' + '+Value.a.a.b.ToDecimalString+' u + '+Value.a.b.a.ToDecimalString+' v + '+Value.a.b.b.ToDecimalString+' uv + '
+Value.a.c.a.ToDecimalString+' v² + '+Value.a.c.b.ToDecimalString+' uv² + '+Value.b.a.a.ToDecimalString+' w + '+Value.b.a.b.ToDecimalString+' uw + '+
Value.b.b.a.ToDecimalString+' vw + '+Value.b.b.b.ToDecimalString+' uvw + '+Value.b.c.a.ToDecimalString+' v²w + '+Value.b.c.b.ToDecimalString+' uv²w';}
end;

        {****** Convert an FP12 to a Hexadecimal String *******}
function _FP12_To_HexString(const value:Fp12Int):string;
begin
Result:=Value.a.ToHexString+'('+value.b.ToHexString+')* w';
{Result:=Value.a.a.a.ToHexString+' + '+Value.a.a.b.ToHexString+' u + '+Value.a.b.a.ToHexString+' v + '+Value.a.b.b.ToHexString+' uv + '
+Value.a.c.a.ToHexString+' v² + '+Value.a.c.b.ToHexString+' uv² + '+Value.b.a.a.ToHexString+' w + '+Value.b.a.b.ToHexString+' uw + '+
Value.b.b.a.ToHexString+' vw + '+Value.b.b.b.ToHexString+' uvw + '+Value.b.c.a.ToHexString+' v²w + '+Value.b.c.b.ToHexString+' uv²w';}
end;

        {****** Convert String to an FP12 (Decimal/Hex)*******}
procedure _FP6_From_Strings(const valueA1,valueB1,valueC1,valueA2,valueB2,valueC2:String;var Result:Fp12Int);
var i:integer;
    s1,s2,Val:string;
    valid:boolean;
begin
Result.a.a.SetFormString(valueA1);
Result.a.b.SetFormString(valueB1);
Result.a.c.SetFormString(valueC1);
Result.b.a.SetFormString(valueA2);
Result.b.b.SetFormString(valueB2);
Result.b.c.SetFormString(valueC2);
end;

{*******************************************************************************}
      ///  Definitions of an FP12 integer operators and functions
{*******************************************************************************}
class operator Fp12Int.Add(const Left, Right: Fp12Int): Fp12Int;
begin
_Add_FP12(Left,Right,Result);
end;

class operator Fp12Int.Subtract(const Left, Right: Fp12Int): Fp12Int;
begin
_Sub_FP12(Left,Right,Result);
end;

class operator Fp12Int.Multiply(const Left, Right: Fp12Int): Fp12Int;
begin
_Mul_FP12(Left,Right,Result);
end;

function Fp12Int.Sqr:Fp12Int;
begin
_Sqr_FP12(Self,Result);
end;

function Fp12Int.Conjugate: Fp12Int;
begin
_Conjugate_FP12(Self,Result);
end;

function Fp12Int.Pow(e:LInt;powerMode:FpPoweringMode=pmNormal):Fp12Int;
begin
_Pow_FP12(Self,e,Result,powerMode);
end;

function Fp12Int.PowToP: Fp12Int;
begin
_Pow_FP12_P(Self,Result);
end;

function Fp12Int.PowToP2: Fp12Int;
begin
_Pow_FP12_P2(Self,Result);
end;

function Fp12Int.PowToP3: Fp12Int;
begin
_Pow_FP12_P3(Self,Result);
end;

function Fp12Int.Inverse: Fp12Int;
begin
_Inv_FP12(Self,Result);
end;

function Fp12Int.IsOne: boolean;
begin
Result:=(a.IsOne)and(b.IsZero);
end;

function Fp12Int.IsZero: boolean;
begin
Result:=(a.IsZero)and(b.IsZero);
end;

class operator Fp12Int.Equal(const Left, Right: Fp12Int): Boolean;
begin
result:=_Equals_FP12(Left,Right);
end;

class operator Fp12Int.NotEqual(const Left, Right: Fp12Int): Boolean;
begin
result:=not _Equals_FP12(Left,Right);
end;

procedure Fp12Int.SetTowerParams(param: PtrTowerParams12);
begin
Tower:=param;
a.SetTowerParams(Tower);
b.SetTowerParams(Tower);
end;

procedure Fp12Int.SetToZero;
begin
a.a.a:=0;a.a.b:=0;
a.b.a:=0;a.b.b:=0;
a.c.a:=0;a.c.b:=0;
b.a.a:=0;b.a.b:=0;
b.b.a:=0;b.b.b:=0;
b.c.a:=0;b.c.b:=0;
end;

procedure Fp12Int.SetFromStrings(a1, b1, c1, a2, b2, c2: String);
begin
_FP6_From_Strings(a1,b1,c1,a2,b2,c2,self);
end;

procedure Fp12Int.SetToOne;
begin
a.a.a:=1;a.a.b:=0;
a.b.a:=0;a.b.b:=0;
a.c.a:=0;a.c.b:=0;
b.a.a:=0;b.a.b:=0;
b.b.a:=0;b.b.b:=0;
b.c.a:=0;b.c.b:=0;
end;

procedure Fp12Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
end;

function Fp12Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.a.BitLength mod 8) * 8;
Setlength(Result,size*12);
Move(a.a.a.Data.i8[0],Result[0],size);
Move(a.a.b.Data.i8[0],Result[size],size);
Move(a.b.a.Data.i8[0],Result[2*size],size);
Move(a.b.b.Data.i8[0],Result[3*size],size);
Move(a.c.a.Data.i8[0],Result[4*size],size);
Move(a.c.b.Data.i8[0],Result[5*size],size);

Move(b.a.a.Data.i8[0],Result[6*Size],size);
Move(b.a.b.Data.i8[0],Result[7*size],size);
Move(b.b.a.Data.i8[0],Result[8*size],size);
Move(b.b.b.Data.i8[0],Result[9*size],size);
Move(b.c.a.Data.i8[0],Result[10*size],size);
Move(b.c.b.Data.i8[0],Result[11*size],size);
end;

function Fp12Int.ToDecimalString: string;
begin
Result:=_FP12_To_DecimalString(Self)
end;

function Fp12Int.ToHexString: string;
begin
Result:=_FP12_To_HexString(Self)
end;

end.
