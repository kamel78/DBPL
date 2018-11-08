unit Fp18Arithmetic;

interface

uses System.SysUtils,VCL.dialogs, LargeIntegers, Fp3Arithmetic,Fp9Arithmetic, GeneralTypes;

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp18. Fp18 is the 18th extention of Fp (Tower extention of order 2 for FP9)
   with respect to the irreductible polynômial W^2-Gamma=0. Elements are in the Form a+b*W (a and b are from Fp9).

********************************************************************************************}

Type



  PFp18Int=^Fp18Int;

  Fp18Int=record
  a,b:Fp9Int;
  Tower:PtrTowerParams18;
  public
    class operator Add(const Left, Right: Fp18Int): Fp18Int;
    class operator Subtract(const Left, Right: Fp18Int): Fp18Int;
    class operator Multiply(const Left, Right: Fp18Int): Fp18Int;
    class operator Equal(const Left, Right: Fp18Int): Boolean;
    class operator NotEqual(const Left, Right: Fp18Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse:Fp18Int;
    function Sqr:Fp18Int;
    function Conjugate:Fp18Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function Pow(e:LInt):Fp18Int;
    function PowToP:Fp18Int;                      /// optimized power computation using Frobenius Map
    function PowToPi(i:integer):Fp18Int;
    procedure SetTowerParams(param:PtrTowerParams18);
    procedure SetToOne;
    procedure SetToZero;
    procedure SetToRandom;
    procedure SetFromStrings(a1,b1,c1,a2,b2,c2:String);
  end;

    procedure _FP18_From_Strings(const valueA1,valueB1,valueC1,valueA2,valueB2,valueC2:String;var Result:Fp18Int);
    function _FP18_To_DecimalString(const value:Fp18Int):string;
    procedure _Inv_FP18(const value:Fp18Int;var Result:Fp18Int);
    procedure _Pow_FP18_P(const value:Fp18Int;var Result:Fp18Int);
    procedure _Pow_FP18_P_i(const value:Fp18Int; power: integer;var Result:Fp18Int);
    procedure _Pow_FP18(const value:Fp18Int;Exponent:LInt;var Result:Fp18Int;Poweringmode :FpPoweringMode=pmNormal);
    procedure _Conjugate_FP18(const Value:Fp18Int;var result:Fp18Int);
    procedure _Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
    procedure _Be_Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
    procedure _Mul_FP18(const Left,Right:Fp18Int; var Result:Fp18Int);
    procedure _Sub_FP18(const Left,Right:Fp18Int; var Result:Fp18Int);
    procedure _Add_FP18(const Left,Right:Fp18Int; var Result:Fp18Int);
    function _Equals_FP18(const Left,Right:Fp18Int):boolean;
    procedure _Compressed_Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
    procedure _DeCompressed_Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
    procedure _Sparse_Mul_FP18(const Left,Right:Fp18Int; var Result:Fp18Int;TwistMode:TTwistModel);

implementation

{*******************************************************************************}
            ///      Procedures for Fp18 Arithmetic
{*******************************************************************************}

        {**********   Add two Fp18 integers *****************}
procedure _Add_FP18(const Left,Right:Fp18Int; var Result:Fp18Int);
begin
Result.Tower:=Left.Tower;
_Add_FP9(Left.a,Right.a,Result.a);
_Add_FP9(Left.b,Right.b,Result.b);
end;

        {**********   Sub two Fp18 integers *****************}
procedure _Sub_FP18(const Left,Right:Fp18Int; var Result:Fp18Int);
begin
Result.Tower:=Left.Tower;
_Sub_FP9(Left.a,Right.a,Result.a);
_Sub_FP9(Left.b,Right.b,Result.b);
end;

        {********** Multiply two Fp18 integers *****************}
procedure _Sparse_Mul_FP18(const Left,Right:Fp18Int; var Result:Fp18Int; TwistMode:TTwistModel);
var t:array[0..4] of Fp9Int;
    Sigmaisiplus1:boolean;
begin
Sigmaisiplus1:=(Left.Tower.Sigma.a.Data.i16[-2]=1)and(LEft.Tower.Sigma.b.Data.i16[-2]=1)and(Left.Tower.Sigma.a.Data.i16[0]=1)and(Left.Tower.Sigma.b.Data.i16[0]=1);
Result.Tower:=Left.Tower;
if TwistMode=twDType then begin
                          _Sparse1_Mul_FP9(Left.a,Right.a,t[0]);
                          _Sparse2_Mul_FP9(Left.b,Right.b,t[1]);
                          end
else begin
     _Sparse2_Mul_FP9(Left.a,Right.a,t[0]);
     _Sparse3_Mul_FP9(Left.b,Right.b,t[1]);
     end;
_Add_FP9(Left.a,Left.b,t[2]);
if TwistMode=twDType then _Sparse_Add_FP9(Right.a,Right.b,t[3])
else _Add_FP9(Right.b,Right.a,t[3]);
if Sigmaisiplus1 then begin
                      /// Works Only if Sigma is 1+i
                      t[4].Tower:=t[1].Tower;
                      t[4].b:=t[1].a;
                      t[4].c:=t[1].b;
                      _Sub_Lint(t[1].c.a,t[1].c.b,t[4].a.a);
                      _Add_Lint(t[1].c.a,t[1].c.b,t[4].a.b);
                      ///
                      end
else _Mul_FP9_By_V(t[1],t[4]);
_Add_FP9(t[0],t[4],Result.a);
_Sparse2_Mul_FP9(t[2],t[3],Result.b);
_Nomod_Sub_FP9(Result.b,t[0],Result.b);
_Sub_FP9(Result.b,t[1],Result.b);
end;

        {********** Multiply two Fp18 integers *****************}
procedure _Mul_FP18(const Left,Right:Fp18Int; var Result:Fp18Int);
var t:array[0..4] of Fp9Int;
begin
Result.Tower:=Left.Tower;
_Mul_FP9(Left.a,Right.a,t[0]);
_Mul_FP9(Left.b,Right.b,t[1]);
_Add_FP9(Left.a,Left.b,t[2]);
_Add_FP9(Right.a,Right.b,t[3]);
_Mul_FP9_By_V(t[1],t[4]);
_Add_FP9(t[0],t[4],Result.a);
_Mul_FP9(t[2],t[3],Result.b);
_Sub_FP9(Result.b,t[0],Result.b);
_Sub_FP9(Result.b,t[1],Result.b);
end;

        {********** Get Square of an Fp18 *****************}
procedure _Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
var  t:array[0..1] of Fp9Int;
      Sigmaisiplus1:boolean;
begin
      {       Using classical Squering    }
Result.Tower:=Value.Tower;
_Sqr_FP9(Value.a,t[0]);
_Sqr_FP9(Value.b,t[1]);
_Mul_FP9(value.a,Value.b,Result.b);
_Mul_FP9_By_V(t[1],Result.a);
_Add_FP9(Result.a,t[0],Result.a);
_Add_FP9(Result.b,Result.b,Result.b);
end;

        {********** Get Square of an FP4 element *****************}
        // Alternative Defined for Beuchat's Algorithm ovec Cyclotomic Group
procedure _Sqr_Fp4(const a0,a1:FP3Int;tower:PtrTowerParams18; var c0,c1:FP3Int);
var t:array[0..1] of FP3Int;
begin
_Sqr_FP3(a0,t[0]);
_Sqr_FP3(a1,t[1]);
_Mul_FP3(t[1],tower^.Sigma,c0);
_Add_FP3(c0,t[0],c0);
_Add_FP3(a0,a1,c1);
_Sqr_FP3(c1,c1);
_Sub_FP3(c1,t[0],c1);
_Sub_FP3(c1,t[1],c1);
end;

     {********** Get Compressed Square of an Fp18 (Karabina Approach)*****************}
     // Value should verify the condition (Value)^(p^6+1)=1 (Cyclotomic)
procedure _Compressed_Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
var  t:array[0..6] of FP3Int;
begin
_Sqr_FP3(value.a.b,t[0]);
_Sqr_FP3(value.b.c,t[1]);
_Add_FP3(value.a.b,value.b.c,t[5]);
_Sqr_FP3(t[5],t[2]);
_Nomod_Add_FP3(t[0],t[1],t[3]);
_Nomod_Sub_FP3(t[2],t[3],t[3]);
_Mod_FP3_FP(t[3],value.Tower.FieldParam.p,t[5]);
_Add_FP3(value.b.a,value.a.c,t[6]);
_Sqr_FP3(t[6],t[3]);
_Sqr_FP3(value.b.a,t[2]);
_Mul_FP3(t[5],value.Tower.Sigma,t[6]);
_Nomod_Add_FP3(t[6],value.b.a,t[5]);
_Nomod_Add_FP3(t[5],t[5],t[5]);
_Add_FP3(t[5],t[6],Result.b.a);
_Mul_FP3(t[1],value.Tower.Sigma,t[4]);
_Nomod_Add_FP3(t[0],t[4],t[5]);
_Nomod_Sub_FP3(t[5],value.a.c,t[6]);
_Add_FP3(t[6],t[6],t[6]);
_Sqr_FP3(value.a.c,t[1]);
_Add_FP3(t[5],t[6],result.a.c);
_Mul_FP3(t[1],value.Tower.Sigma,t[4]);
_Add_FP3(t[4],t[2],t[5]);
_Nomod_Sub_FP3(t[5],value.a.b,t[6]);
_Add_FP3(t[6],t[6],t[6]);
_Add_FP3(t[5],t[6],result.a.b);
_Nomod_Add_FP3(t[1],t[2],t[0]);
_Nomod_Sub_FP3(t[3],t[0],t[5]);
_Nomod_Add_FP3(t[5],value.b.c,t[6]);
_Nomod_Add_FP3(t[6],t[6],t[6]);
_Add_FP3(t[6],t[5],Result.b.c);
Result.b.b.a:=0;
Result.b.b.b:=0;
Result.a.a.a:=0;
Result.a.a.b:=0;
end;

    {********** Get Decompressed Square of an Fp18 (Karabina Approach)*****************}
procedure _DeCompressed_Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
var t:array[0..6] of FP3Int;
begin
if _Is_FP3_Null(value.a.b) then begin
                                _Mul_FP3(value.a.b,value.b.c,t[1]);
                                _Add_FP3(t[1],t[1],t[1]);
                                _Inv_FP3(value.a.c,t[2]);
                                _Mul_FP3(t[1],t[2],result.b.b);
                                _Sqr_FP3(t[1],t[2]);
                                _Add_FP3(t[2],t[2],t[2]);
                                _Mul_FP3(value.a.c,value.a.b,t[3]);
                                _Add_FP3(t[3],t[3],t[4]);
                                _Add_FP3(t[3],t[4],t[4]);
                                _Sub_FP3(t[3],t[4],t[0]);
                                _Mul_FP3(t[0],value.Tower.Sigma,Result.a.a);
                                _Inc_LInt(Result.a.a.a,1);
                                end
else begin
     _Sqr_FP3(value.b.c,t[5]);
     _Mul_FP3(t[5],value.Tower.Sigma,t[5]);
     _Sqr_FP3(value.a.b,t[4]);
     _Add_FP3(t[4],t[4],t[6]);
     _Add_FP3(t[4],t[6],t[4]);
     _Add_FP3(value.a.c,value.a.c,t[3]);
     _Add_FP3(value.b.a,value.b.a,t[2]);
     _Add_FP3(t[2],t[2],t[2]);
     _Inv_FP3(t[2],t[2]);
     _Add_FP3(t[4],t[5],t[1]);
     _Sub_FP3(t[1],t[3],t[1]);
     _Mul_FP3(t[1],t[2],Result.b.b);
     _Sqr_FP3(Result.b.b,t[1]);
     _Add_FP3(t[1],t[1],t[1]);
     _Mul_FP3(value.b.a,value.b.c,t[2]);
     _Mul_FP3(value.a.c,value.a.b,t[3]);
     _Add_FP3(t[3],t[3],t[4]);
     _Add_FP3(t[3],t[4],t[3]);
     _Add_FP3(t[1],t[2],t[0]);
     _Sub_FP3(t[0],t[3],t[0]);
     _Mul_FP3(t[0],value.Tower.Sigma,Result.a.a);
     _Inc_LInt(Result.a.a.a,1);
     end;
Result.a.b:=Value.a.b;
Result.a.c:=Value.a.c;
Result.b.a:=Value.b.a;
Result.a.b:=Value.a.b;
Result.b.c:=Value.b.c;
end;

        {********** Get Square of an Fp18 (Beuchat Approach)*****************}
          // Value should verify the condition (Value)^(p^6+1)=1
procedure _Be_Sqr_FP18(const value: Fp18Int; var Result: Fp18Int);
var  t00,t11,t10,t01,t02,t12,aux,tmp: FP3Int;
begin
      {       Using Beuchat FP4 technique    }
Result.Tower:=Value.Tower;
_Sqr_Fp4(value.a.a,value.b.b,value.Tower,t00,t11);
_Sqr_Fp4(value.b.a,value.a.c,value.Tower,t01,t12);
_Sqr_Fp4(value.a.b,value.b.c,value.Tower,t02,aux);
_Mul_FP3(aux,value.Tower^.Sigma,t10);
_Add_FP3(t00,t00,tmp);
_Add_FP3(t00,tmp,tmp);
_Sub_FP3(tmp,value.a.a,tmp);
_Sub_FP3(tmp,value.a.a,Result.a.a);
_Add_FP3(t01,t01,tmp);
_Add_FP3(t01,tmp,tmp);
_Sub_FP3(tmp,value.a.b,tmp);
_Sub_FP3(tmp,value.a.b,Result.a.b);
_Add_FP3(t02,t02,tmp);
_Add_FP3(t02,tmp,tmp);
_Sub_FP3(tmp,value.a.c,tmp);
_Sub_FP3(tmp,value.a.c,Result.a.c);
_Add_FP3(t10,t10,tmp);
_Add_FP3(t10,tmp,tmp);
_Add_FP3(tmp,value.b.a,tmp);
_Add_FP3(tmp,value.b.a,Result.b.a);
_Add_FP3(t11,t11,tmp);
_Add_FP3(t11,tmp,tmp);
_Add_FP3(tmp,value.b.b,tmp);
_Add_FP3(tmp,value.b.b,Result.b.b);
_Add_FP3(t12,t12,tmp);
_Add_FP3(t12,tmp,tmp);
_Add_FP3(tmp,value.b.c,tmp);
_Add_FP3(tmp,value.b.c,Result.b.c);
end;

        {**********   Conjugate an Fp18 integer *************}
procedure _Conjugate_FP18(const Value:Fp18Int;var result:Fp18Int);
begin
Result.Tower:=Value.Tower;
Result.a:=Value.a;
_Neg_FP9(Value.b,Result.b);
end;

        {********** Raise an Fp18 to a LInt power ************}
procedure _Pow_FP18(const value:Fp18Int;Exponent:LInt;var Result:Fp18Int;Poweringmode :FpPoweringMode=pmNormal);
var i:word;
    zi,tmp:Fp18Int;
    nafexpo,t1,t2:LIntArrayForm;
    SizeNaf,Sizebin:Integer;
begin
Result.Tower:=Value.Tower;
case Poweringmode of
pmNormal:begin
              //  Raise an Fp18 to a LInt power (Classical Square and Multiply approach)************
         Result:=Value;
         for i:=Exponent.BitLength-2 downto 0 do begin
                                                 _Sqr_FP18(Result,Result);
                                                 if _Is_BitSet_At(Exponent,i) then _Mul_FP18(Result,value,Result);
                                                 end;
         end;
pmBeuchat:begin
                  //  Raise an Fp18 to a LInt power (Beucaht FP4 Approach)************
                  //  Value should verify the condition (Value)^(p^6+1)=1
         Result:=Value;
         for i:=Exponent.BitLength-2 downto 0 do begin
                                                 _Be_Sqr_FP18(Result,Result);
                                                 if _Is_BitSet_At(Exponent,i) then _Mul_FP18(Result,value,Result);
                                                 end;
          end;
pmKarbina:begin
                  //  Raise an Fp18 to a LInt power (Karabina Approach)************
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
                                     if NafExpo[0]=-1 then _Conjugate_FP18(Result,Result);
                                     end
          else begin
               Result.SetToOne;
               Result.SetTowerParams(Value.Tower);
               end;
          for i:=1 to Length(NafExpo)-1 do begin
                                           _Compressed_Sqr_FP18(tmp,tmp);
                                           if nafexpo[i]=1 then begin
                                                                _DeCompressed_Sqr_FP18(tmp,zi);
                                                                _Mul_FP18(zi,Result,Result);
                                                                end
                                           else if nafexpo[i]=-1 then begin
                                                                      _DeCompressed_Sqr_FP18(tmp,zi);
                                                                      _Conjugate_FP18(zi,zi);
                                                                      _Mul_FP18(zi,Result,Result);
                                                                      end
                                           end;
          end;
end;
if _IsNeg(Exponent) then _Conjugate_FP18(Result,Result);
end;

        {********** Raise an Fp18 to " p " power ************}
        /// use Precomputed Frobenius Constants
procedure _Pow_FP18_P(const value:Fp18Int;var Result:Fp18Int);
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.a.a:=Value.a.a.toPowerP;
_Mul_FP3(Value.a.b.toPowerP,Value.Tower^.FrobeniusP_Const[1],Result.a.b);
_Mul_FP3(Value.a.c.toPowerP,Value.Tower^.FrobeniusP_Const[3],Result.a.c);
_Mul_FP3(Value.b.a.toPowerP,Value.Tower^.FrobeniusP_Const[0],Result.b.a);
_Mul_FP3(Value.b.b.toPowerP,Value.Tower^.FrobeniusP_Const[2],Result.b.b);
_Mul_FP3(Value.b.c.toPowerP,Value.Tower^.FrobeniusP_Const[4],Result.b.c);
end;

        {********** Raise an Fp18 to " p^i " power ************}
        /// use Precomputed Frobenius Constant
procedure _Pow_FP18_P_i(const value:Fp18Int; power: integer;var Result:Fp18Int);
var i:integer;
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result := Value;
for i := 0 to power - 1 do _Pow_FP18_P(Result, Result);
end;



        {********** Get inverse of an Fp18 *****************}
procedure _Inv_FP18(const value:Fp18Int;var Result:Fp18Int);
var t:array[0..2]of Fp9Int;
begin
Result.Tower:=Value.Tower;
_Sqr_FP9(value.a,t[0]);
_Sqr_FP9(value.b,t[1]);
_Mul_FP9_By_V(t[1],t[2]);
_Sub_FP9(t[0],t[2],t[0]);
_Inv_FP9(t[0],t[1]);
_Mul_FP9(Value.a,t[1],Result.a);
_Mul_FP9(Value.b,t[1],Result.b);
_Neg_FP9(Result.b,Result.b);
end;

        {**********   Compare two Fp18 integers *************}
function _Equals_FP18(const Left,Right:Fp18Int):boolean;
begin
result:=_Equals_FP9(Left.a,Right.a) and _Equals_FP9(Left.b,Right.b);
end;

        {****** Convert an Fp18 to a Decimal String **********}
function _FP18_To_DecimalString(const value:Fp18Int):string;
begin
Result:=Value.a.a.a.ToDecimalString+' + '+Value.a.a.b.ToDecimalString+' u + '+Value.a.a.c.ToDecimalString+' u^2 + '+
Value.a.b.a.ToDecimalString+' v + '+Value.a.b.b.ToDecimalString+' uv + '+Value.a.b.c.ToDecimalString+' u^2v + '
+Value.a.c.a.ToDecimalString+' v^2 + '+Value.a.c.b.ToDecimalString+' uv^2 + '+Value.a.c.c.ToDecimalString+' u^2v^2 + '
+Value.b.a.a.ToDecimalString+' w + '+Value.b.a.b.ToDecimalString+' uw + '+Value.b.a.c.ToDecimalString+' u^2w + '+
Value.b.b.a.ToDecimalString+' vw + '+Value.b.b.b.ToDecimalString+' uvw + '+Value.b.b.c.ToDecimalString+' u^2vw + '
+Value.b.c.a.ToDecimalString+' v^2w + '+Value.b.c.b.ToDecimalString+' uv^2w + '+Value.b.c.c.ToDecimalString+' u^2v^2w';
end;

        {****** Convert an Fp18 to a Hexadecimal String *******}
function _FP18_To_HexString(const value:Fp18Int):string;
begin
Result:=Value.a.a.a.ToHexString+' + '+Value.a.a.b.ToHexString+' u + '+Value.a.a.c.ToHexString+' u^2 + '+
Value.a.b.a.ToHexString+' v + '+Value.a.b.b.ToHexString+' uv + '+Value.a.b.c.ToHexString+' u^2v + '
+Value.a.c.a.ToHexString+' v^2 + '+Value.a.c.b.ToHexString+' uv^2 + '+Value.a.c.c.ToHexString+' u^2v^2 + '
+Value.b.a.a.ToHexString+' w + '+Value.b.a.b.ToHexString+' uw + '+Value.b.a.c.ToHexString+' u^2w + '+
Value.b.b.a.ToHexString+' vw + '+Value.b.b.b.ToHexString+' uvw + '+Value.b.b.c.ToHexString+' u^2vw + '
+Value.b.c.a.ToHexString+' v^2w + '+Value.b.c.b.ToHexString+' uv^2w + '+Value.b.c.c.ToHexString+' u^2v^2w';
end;

        {****** Convert String to an Fp18 (Decimal/Hex)*******}
procedure _FP18_From_Strings(const valueA1,valueB1,valueC1,valueA2,valueB2,valueC2:String;var Result:Fp18Int);
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
      ///  Definitions of an Fp18 integer operators and functions
{*******************************************************************************}
class operator Fp18Int.Add(const Left, Right: Fp18Int): Fp18Int;
begin
_Add_FP18(Left,Right,Result);
end;

class operator Fp18Int.Subtract(const Left, Right: Fp18Int): Fp18Int;
begin
_Sub_FP18(Left,Right,Result);
end;

class operator Fp18Int.Multiply(const Left, Right: Fp18Int): Fp18Int;
begin
_Mul_FP18(Left,Right,Result);
end;

function Fp18Int.Sqr:Fp18Int;
begin
_Sqr_FP18(Self,Result);
end;

function Fp18Int.Conjugate: Fp18Int;
begin
_Conjugate_FP18(Self,Result);
end;

function Fp18Int.Pow(e: LInt): Fp18Int;
begin
_Pow_FP18(Self,e,Result);
end;

function Fp18Int.PowToP: Fp18Int;
begin
_Pow_FP18_P(Self,Result);
end;

function Fp18Int.PowToPi(i:integer): Fp18Int;
begin
_Pow_FP18_P_i(Self,i,result);
end;

function Fp18Int.Inverse: Fp18Int;
begin
_Inv_FP18(Self,Result);
end;

function Fp18Int.IsOne: boolean;
begin
Result:=(a.IsOne)and(b.IsZero);
end;

function Fp18Int.IsZero: boolean;
begin
Result:=(a.IsZero)and(b.IsZero);
end;

class operator Fp18Int.Equal(const Left, Right: Fp18Int): Boolean;
begin
result:=_Equals_FP18(Left,Right);
end;

class operator Fp18Int.NotEqual(const Left, Right: Fp18Int): Boolean;
begin
result:=not _Equals_FP18(Left,Right);
end;

procedure Fp18Int.SetTowerParams(param: PtrTowerParams18);
begin
Tower:=param;
a.SetTowerParams(Tower);
b.SetTowerParams(Tower);
end;

procedure Fp18Int.SetToZero;
begin
a.SetToZero;
b.SetToZero;
end;

procedure Fp18Int.SetFromStrings(a1, b1, c1, a2, b2, c2: String);
begin
_FP18_From_Strings(a1,b1,c1,a2,b2,c2,self);
end;

procedure Fp18Int.SetToOne;
begin
a.SetToOne;
b.SetToZero;
end;

procedure Fp18Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
end;

function Fp18Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.a.BitLength mod 8) * 8;
Setlength(Result,size*18);
Move(a.a.a.Data.i8[0],Result[0],size);
Move(a.a.b.Data.i8[0],Result[size],size);
Move(a.a.c.Data.i8[0],Result[2*size],size);
Move(a.b.a.Data.i8[0],Result[3*size],size);
Move(a.b.b.Data.i8[0],Result[4*size],size);
Move(a.b.c.Data.i8[0],Result[5*size],size);
Move(a.c.a.Data.i8[0],Result[6*size],size);
Move(a.c.b.Data.i8[0],Result[7*size],size);
Move(a.c.c.Data.i8[0],Result[8*size],size);

Move(b.a.a.Data.i8[0],Result[9*size],size);
Move(b.a.b.Data.i8[0],Result[10*size],size);
Move(b.a.c.Data.i8[0],Result[11*size],size);
Move(b.b.a.Data.i8[0],Result[12*size],size);
Move(b.b.b.Data.i8[0],Result[13*size],size);
Move(b.b.c.Data.i8[0],Result[14*size],size);
Move(b.c.a.Data.i8[0],Result[15*size],size);
Move(b.c.b.Data.i8[0],Result[16*size],size);
Move(b.c.c.Data.i8[0],Result[17*size],size);
end;

function Fp18Int.ToDecimalString: string;
begin
Result:=_FP18_To_DecimalString(Self)
end;

function Fp18Int.ToHexString: string;
begin
Result:=_FP18_To_HexString(Self)
end;

end.
