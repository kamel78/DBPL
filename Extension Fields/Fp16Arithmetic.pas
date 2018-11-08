unit Fp16Arithmetic;

interface

uses System.SysUtils,VCL.dialogs, LargeIntegers, Fp4Arithmetic,Fp8Arithmetic, GeneralTypes;

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp16. Fp16 is the 16th extention of Fp (Tower extention of order 2 for FP8)
   with respect to the irreductible polynômial W^2-Gamma=0. Elements are in the Form a+b*W (a and b are from Fp8).

********************************************************************************************}

Type



  PFp16Int=^Fp16Int;

  Fp16Int=record
  a,b:Fp8Int;
  Tower:PtrTowerParams24;
  public
    class operator Add(const Left, Right: Fp16Int): Fp16Int;
    class operator Subtract(const Left, Right: Fp16Int): Fp16Int;
    class operator Multiply(const Left, Right: Fp16Int): Fp16Int;
    class operator Equal(const Left, Right: Fp16Int): Boolean;
    class operator NotEqual(const Left, Right: Fp16Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse:Fp16Int;
    function Sqr:Fp16Int;
    function Conjugate:Fp16Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function Pow(e:LInt):Fp16Int;
    function PowToP:Fp16Int;                      /// optimized power computation using Frobenius Map
    function PowToPi(i:integer):Fp16Int;
    procedure SetTowerParams(param:PtrTowerParams24);
    procedure SetToOne;
    procedure SetToZero;
    procedure SetToRandom;
    procedure SetFromStrings(a1,b1,c1,a2,b2,c2:String);
  end;

    procedure _FP16_From_Strings(const valueA1,valueB1,valueC1,valueA2,valueB2,valueC2:String;var Result:Fp16Int);
    function _FP16_To_DecimalString(const value:Fp16Int):string;
    procedure _Inv_FP16(const value:Fp16Int;var Result:Fp16Int);
    procedure _Pow_FP16_P(const value:Fp16Int;var Result:Fp16Int);
    procedure _Pow_FP16_P_i(const value:Fp16Int; power: integer;var Result:Fp16Int);
    procedure _Pow_FP16(const value:Fp16Int;Exponent:LInt;var Result:Fp16Int;Poweringmode :FpPoweringMode=pmNormal);
    procedure _Conjugate_FP16(const Value:Fp16Int;var result:Fp16Int);
    procedure _Sqr_FP16(const value: Fp16Int; var Result: Fp16Int);
    procedure _Mul_FP16(const Left,Right:Fp16Int; var Result:Fp16Int);
    procedure _Sub_FP16(const Left,Right:Fp16Int; var Result:Fp16Int);
    procedure _Add_FP16(const Left,Right:Fp16Int; var Result:Fp16Int);
    function _Equals_FP16(const Left,Right:Fp16Int):boolean;

implementation

{*******************************************************************************}
            ///      Procedures for Fp16 Arithmetic
{*******************************************************************************}

        {**********   Add two Fp16 integers *****************}
procedure _Add_FP16(const Left,Right:Fp16Int; var Result:Fp16Int);
begin
Result.Tower:=Left.Tower;
_Add_FP8(Left.a,Right.a,Result.a);
_Add_FP8(Left.b,Right.b,Result.b);
end;

        {**********   Sub two Fp16 integers *****************}
procedure _Sub_FP16(const Left,Right:Fp16Int; var Result:Fp16Int);
begin
Result.Tower:=Left.Tower;
_Sub_FP8(Left.a,Right.a,Result.a);
_Sub_FP8(Left.b,Right.b,Result.b);
end;

        {********** Multiply two Fp16 integers *****************}
{procedure _Sparse_Mul_FP16(const Left,Right:Fp16Int; var Result:Fp16Int; TwistMode:TTwistModel);
var t:array[0..4] of Fp8Int;
    Sigmaisiplus1:boolean;
begin
Sigmaisiplus1:=(Left.Tower.Sigma.a.Data.i16[-2]=1)and(LEft.Tower.Sigma.b.Data.i16[-2]=1)and(Left.Tower.Sigma.a.Data.i16[0]=1)and(Left.Tower.Sigma.b.Data.i16[0]=1);
Result.Tower:=Left.Tower;
if TwistMode=twDType then begin
                          _Sparse1_Mul_FP8(Left.a,Right.a,t[0]);
                          _Sparse2_Mul_FP8(Left.b,Right.b,t[1]);
                          end
else begin
     _Sparse2_Mul_FP8(Left.a,Right.a,t[0]);
     _Sparse3_Mul_FP8(Left.b,Right.b,t[1]);
     end;
_Add_FP8(Left.a,Left.b,t[2]);
if TwistMode=twDType then _Sparse_Add_FP8(Right.a,Right.b,t[3])
else _Add_FP8(Right.b,Right.a,t[3]);
if Sigmaisiplus1 then begin
                      /// Works Only if Sigma is 1+i
                      t[4].Tower:=t[1].Tower;
                      t[4].b:=t[1].a;
                      t[4].c:=t[1].b;
                      _Sub_Lint(t[1].c.a,t[1].c.b,t[4].a.a);
                      _Add_Lint(t[1].c.a,t[1].c.b,t[4].a.b);
                      ///
                      end
else _Mul_FP8_By_V(t[1],t[4]);
_Add_FP8(t[0],t[4],Result.a);
_Sparse2_Mul_FP8(t[2],t[3],Result.b);
_Nomod_Sub_FP8(Result.b,t[0],Result.b);
_Sub_FP8(Result.b,t[1],Result.b);
end;    }

        {********** Multiply two Fp16 integers *****************}
procedure _Mul_FP16(const Left,Right:Fp16Int; var Result:Fp16Int);
var t:array[0..4] of Fp8Int;
begin
Result.Tower:=Left.Tower;
_Mul_FP8(Left.a,Right.a,t[0]);
_Mul_FP8(Left.b,Right.b,t[1]);
_Add_FP8(Left.a,Left.b,t[2]);
_Add_FP8(Right.a,Right.b,t[3]);
_Mul_FP8_By_w(t[1],t[4]);
_Add_FP8(t[0],t[4],Result.a);
_Mul_FP8(t[2],t[3],Result.b);
_Sub_FP8(Result.b,t[0],Result.b);
_Sub_FP8(Result.b,t[1],Result.b);
end;

        {********** Get Square of an Fp16 *****************}
procedure _Sqr_FP16(const value: Fp16Int; var Result: Fp16Int);
var  t:array[0..1] of Fp8Int;
      Sigmaisiplus1:boolean;
begin
      {       Using classical Squering    }
Result.Tower:=Value.Tower;
_Sqr_FP8(Value.a,t[0]);
_Sqr_FP8(Value.b,t[1]);
_Mul_FP8(value.a,Value.b,Result.b);
_Mul_FP8_By_w(t[1],Result.a);
_Add_FP8(Result.a,t[0],Result.a);
_Add_FP8(Result.b,Result.b,Result.b);
end;

        {**********   Conjugate an Fp16 integer *************}
procedure _Conjugate_FP16(const Value:Fp16Int;var result:Fp16Int);
begin
Result.Tower:=Value.Tower;
Result.a:=Value.a;
_Neg_FP8(Value.b,Result.b);
end;

        {********** Raise an Fp16 to a LInt power ************}
procedure _Pow_FP16(const value:Fp16Int;Exponent:LInt;var Result:Fp16Int;Poweringmode :FpPoweringMode=pmNormal);
var i:word;
    tmp,tmp1:Fp16Int;
    nafexpo,t1,t2:LIntArrayForm;
    SizeNaf,Sizebin:Integer;
begin
//  Raise an Fp16 to a LInt power
//  Value should verify the condition (Value)^(p^8+1)=1

Result.Tower:=Value.Tower;
tmp:=value;
t1:=LIntToNAF(Exponent.Absolute);
t2:=LIntToIntArray(Exponent.Absolute);
SizeNaf:=0;
SizeBin:=0;
for i:=0 to length(t1)-1 do if t1[i]<>0 then inc(SizeNaf);
for i:=0 to length(t2)-1 do if t2[i]<>0 then inc(SizeBin);
if SizeNaf>SizeBin then NafExpo:=t2
else NafExpo:=t1;
if NafExpo[0]<>0  then begin
                       Result:=Value;
                       if NafExpo[0]=-1 then _Conjugate_FP16(Result,Result);
                       end
else begin
     Result.SetToOne;
     Result.SetTowerParams(Value.Tower);
     end;
for i:=1 to Length(NafExpo)-1 do begin
                                 _Sqr_FP16(tmp,tmp);
                                 if nafexpo[i]=1 then _Mul_FP16(tmp,Result,Result)
                                 else if nafexpo[i]=-1 then begin
                                                            _Conjugate_FP16(tmp,tmp1);
                                                            _Mul_FP16(tmp1,Result,Result);
                                                            end
                                 end;
if _IsNeg(Exponent) then _Conjugate_FP16(Result,Result);
end;

        {********** Raise an Fp16 to " p " power ************}
        /// use Precomputed Frobenius Constants
procedure _Pow_FP16_P(const value:Fp16Int;var Result:Fp16Int);
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.a.a:=Value.a.a.toPowerP;
_Mul_FP4(Value.a.b.toPowerP,Value.Tower^.FrobeniusP_Const_Fp4[1],Result.a.b);
_Mul_FP4(Value.b.a.toPowerP,Value.Tower^.FrobeniusP_Const_Fp4[0],Result.b.a);
_Mul_FP4(Value.b.b.toPowerP,Value.Tower^.FrobeniusP_Const_Fp4[2],Result.b.b);
end;

        {********** Raise an Fp16 to " p^i " power ************}
        /// use Precomputed Frobenius Constant
procedure _Pow_FP16_P_i(const value:Fp16Int; power: integer;var Result:Fp16Int);
var i:integer;
begin
Result.Tower:=Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result := Value;
for i := 0 to power - 1 do _Pow_FP16_P(Result, Result);
end;



        {********** Get inverse of an Fp16 *****************}
procedure _Inv_FP16(const value:Fp16Int;var Result:Fp16Int);
var t:array[0..2]of Fp8Int;
begin
Result.Tower:=Value.Tower;
_Sqr_FP8(value.a,t[0]);
_Sqr_FP8(value.b,t[1]);
_Mul_FP8_By_w(t[1],t[2]);
_Sub_FP8(t[0],t[2],t[0]);
_Inv_FP8(t[0],t[1]);
_Mul_FP8(Value.a,t[1],Result.a);
_Mul_FP8(Value.b,t[1],Result.b);
_Neg_FP8(Result.b,Result.b);
end;

        {**********   Compare two Fp16 integers *************}
function _Equals_FP16(const Left,Right:Fp16Int):boolean;
begin
result:=_Equals_FP8(Left.a,Right.a) and _Equals_FP8(Left.b,Right.b);
end;

        {****** Convert an Fp16 to a Decimal String **********}
function _FP16_To_DecimalString(const value:Fp16Int):string;
begin
Result:=Value.a.ToDecimalString+' + ('+Value.b.ToDecimalString+')*z';
end;

        {****** Convert an Fp16 to a Hexadecimal String *******}
function _FP16_To_HexString(const value:Fp16Int):string;
begin
Result:=Value.a.ToHexString+' + ('+Value.b.ToHexString+')*z';
end;

        {****** Convert String to an Fp16 (Decimal/Hex)*******}
procedure _FP16_From_Strings(const valueA1,valueB1,valueC1,valueA2,valueB2,valueC2:String;var Result:Fp16Int);
var i:integer;
    s1,s2,Val:string;
    valid:boolean;
begin
Result.a.a.SetFromString(valueA1);
Result.a.b.SetFromString(valueB1);
Result.b.a.SetFromString(valueA2);
Result.b.b.SetFromString(valueB2);
end;

{*******************************************************************************}
      ///  Definitions of an Fp16 integer operators and functions
{*******************************************************************************}
class operator Fp16Int.Add(const Left, Right: Fp16Int): Fp16Int;
begin
_Add_FP16(Left,Right,Result);
end;

class operator Fp16Int.Subtract(const Left, Right: Fp16Int): Fp16Int;
begin
_Sub_FP16(Left,Right,Result);
end;

class operator Fp16Int.Multiply(const Left, Right: Fp16Int): Fp16Int;
begin
_Mul_FP16(Left,Right,Result);
end;

function Fp16Int.Sqr:Fp16Int;
begin
_Sqr_FP16(Self,Result);
end;

function Fp16Int.Conjugate: Fp16Int;
begin
_Conjugate_FP16(Self,Result);
end;

function Fp16Int.Pow(e: LInt): Fp16Int;
begin
_Pow_FP16(Self,e,Result);
end;

function Fp16Int.PowToP: Fp16Int;
begin
_Pow_FP16_P(Self,Result);
end;

function Fp16Int.PowToPi(i:integer): Fp16Int;
begin
_Pow_FP16_P_i(Self,i,result);
end;

function Fp16Int.Inverse: Fp16Int;
begin
_Inv_FP16(Self,Result);
end;

function Fp16Int.IsOne: boolean;
begin
Result:=(a.IsOne)and(b.IsZero);
end;

function Fp16Int.IsZero: boolean;
begin
Result:=(a.IsZero)and(b.IsZero);
end;

class operator Fp16Int.Equal(const Left, Right: Fp16Int): Boolean;
begin
result:=_Equals_FP16(Left,Right);
end;

class operator Fp16Int.NotEqual(const Left, Right: Fp16Int): Boolean;
begin
result:=not _Equals_FP16(Left,Right);
end;

procedure Fp16Int.SetTowerParams(param: PtrTowerParams24);
begin
Tower:=param;
a.SetTowerParams(Tower);
b.SetTowerParams(Tower);
end;

procedure Fp16Int.SetToZero;
begin
a.SetToZero;
b.SetToZero;
end;

procedure Fp16Int.SetFromStrings(a1, b1, c1, a2, b2, c2: String);
begin
_FP16_From_Strings(a1,b1,c1,a2,b2,c2,self);
end;

procedure Fp16Int.SetToOne;
begin
a.SetToOne;
b.SetToZero;
end;

procedure Fp16Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
end;

function Fp16Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.a.a.BitLength mod 8) * 8;
Setlength(Result,size*16);
Move(a.a.a.a.Data.i8[0],Result[0],size);
Move(a.a.a.b.Data.i8[0],Result[size],size);
Move(a.a.b.a.Data.i8[0],Result[2*size],size);
Move(a.a.b.b.Data.i8[0],Result[3*size],size);
Move(a.b.a.a.Data.i8[0],Result[4*size],size);
Move(a.b.a.b.Data.i8[0],Result[5*size],size);
Move(a.b.b.a.Data.i8[0],Result[6*size],size);
Move(a.b.b.b.Data.i8[0],Result[7*size],size);

Move(b.a.a.a.Data.i8[0],Result[8*size],size);
Move(b.a.a.b.Data.i8[0],Result[9*size],size);
Move(b.a.b.a.Data.i8[0],Result[10*size],size);
Move(b.a.b.b.Data.i8[0],Result[11*size],size);
Move(b.b.a.a.Data.i8[0],Result[12*size],size);
Move(b.b.a.b.Data.i8[0],Result[13*size],size);
Move(b.b.b.a.Data.i8[0],Result[14*size],size);
Move(b.b.b.b.Data.i8[0],Result[15*size],size);
end;

function Fp16Int.ToDecimalString: string;
begin
Result:=_FP16_To_DecimalString(Self)
end;

function Fp16Int.ToHexString: string;
begin
Result:=_FP16_To_HexString(Self)
end;

end.
