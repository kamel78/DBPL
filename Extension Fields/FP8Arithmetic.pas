unit FP8Arithmetic;

interface

uses System.SysUtils,VCL.dialogs, LargeIntegers, Fp2Arithmetic,Fp4Arithmetic,Fp6Arithmetic;

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp8. Fp8 is  the quadratic extention of Fp4 (Tower extention of order 2 for FP4)
   with respect to the polynomial W^2-Gamma=0 (W^2=Gamma). Elements are in the form a+b*W (a and b are from Fp4).

********************************************************************************************}

type

  PtrTowerParams24=^TowerParams24;
  TowerParams24=record
              FieldParam:PtrFieldParams;                // Parametres of the base Field
              Tower8paprms:PtrTowerParams12;            // Parametres of the Tower Fp2<->Fp4
              Sigma:Fp2Int;                             // Parameter of Fp2->Fp4 Field
              Gamma:Fp4Int;                             // Parameter of Fp4->Fp8,Fp12 Field
              Gamma2:Fp6Int;                             // Parameter of Fp6->Fp36 Field
              pmod8:Lint;                               // For BLS24
              FrobeniusPi3xSigmaSqr{Z^(3p-9)^2}: Fp2Int;  // Frobenius Constants for Computation of f^p, f^p^2, f^p^3 ...f^p^7  for BLS24 Curves
              FrobeniusP_Const_Fp4:array[0..4]of Fp4Int;         /// For BLS24 Curves
              FrobeniusP_Const_Fp6:array[0..4]of Fp6Int;         /// For KSS36 Curves
              end;

  Fp8Int=record
  a,b:Fp4Int;
  Tower:PtrTowerParams24;
  public
    class operator Add(const Left, Right: Fp8Int): Fp8Int;
    class operator Subtract(const Left, Right: Fp8Int): Fp8Int;
    class operator Multiply(const Left, Right: Fp8Int): Fp8Int;
    class operator Multiply(const Left:Fp4Int; Right: Fp8Int): Fp8Int;
    class operator Negative(const Value: Fp8Int): Fp8Int;
    class operator Equal(const Left, Right: Fp8Int): Boolean;
    class operator NotEqual(const Left, Right: Fp8Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse:Fp8Int;
    function Pow(Exponent: LInt): Fp8Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function Sqr:Fp8Int;
    function MultiplyByGamma:Fp8Int;
    procedure SetToRandom;
    procedure SetToOne;
    procedure SetToZero;
    procedure SetTowerParams(TParam:PtrTowerParams24);
    procedure SetFromStrings(a,b,c,d:String);
  end;

  procedure _Add_FP8(const Left,Right:Fp8Int; var Result:Fp8Int);
  procedure _Sub_FP8(Const Left,Right:Fp8Int;var Result: Fp8Int);
  procedure _Mul_FP8(const Left,Right:Fp8Int;var Result:Fp8Int);
  procedure _Mul_FP4_FP8(const Left: Fp4Int; Right: Fp8Int; var Result: Fp8Int);
  procedure _Mul_FP8_By_W(const value: Fp8Int; var Result: Fp8Int);
  procedure _Sqr_FP8(const value: Fp8Int; var Result: Fp8Int);
  function  _Equals_FP8(const Left,Right:Fp8Int):boolean;
  procedure _Inv_FP8(const value:Fp8Int;var Result:Fp8Int);
  procedure _Neg_FP8(const Value: Fp8Int; var Result: Fp8Int);
  procedure _Conjugate_FP8(const Source: Fp8Int; var Result: Fp8Int);
  procedure _Pow_FP8(const value:Fp8Int;Exponent:LInt;var Result:Fp8Int);
  function  _FP8_To_DecimalString(const value:Fp8Int):string;
  function  _FP8_To_HexString(const value:Fp8Int):string;
  procedure _FP8_From_Strings(const valueA,valueB,valueC,valueD:String;var Result:Fp8Int);
  procedure _Mul_FP2_FP8(const Left: Fp2Int; Right: Fp8Int; var Result: Fp8Int);
  procedure _Sparse_Mul_FP8(const Left,Right:Fp8Int;var Result:Fp8Int);


implementation

{*******************************************************************************}
            //      Procedures for FP8 Arithmetic
{*******************************************************************************}

        {**********   Add two FP8 integers *****************}
procedure _Add_FP8(const Left,Right:Fp8Int; var Result:Fp8Int);
begin
Result.Tower:=Left.Tower;
_Add_FP4(Left.a,Right.a,Result.a);
_Add_FP4(Left.b,Right.b,Result.b);
end;
        {**********   Substract two FP8 integers ***********}
procedure _Sub_FP8(Const Left,Right:Fp8Int;var Result: Fp8Int);
begin
Result.Tower:=Left.Tower;
_Sub_FP4(Left.a,Right.a,Result.a);
_Sub_FP4(Left.b,Right.b,Result.b);
end;

        {********** Multiply two FP8 integers *****************}
procedure _Sparse_Mul_FP8(const Left,Right:Fp8Int;var Result:Fp8Int);
          /// Right.b is always 0
begin
Result.Tower:=Left.Tower;
_Mul_FP4(Left.a,Right.a,Result.a);
_Mul_FP4(Left.b,Right.a,Result.b);
end;

        {********** Multiply two FP8 integers *****************}
procedure _Mul_FP8(const Left,Right:Fp8Int;var Result:Fp8Int);
var tmp:array[0..4] of Fp4Int;
begin
        { Using Karatsuba Multiplication  }
Result.Tower:=Left.Tower;
_Mul_FP4(Left.a,Right.a,tmp[0]);
_Mul_FP4(Left.b,Right.b,tmp[1]);
_Add_FP4(Left.a,Left.b,tmp[3]);
_Add_FP4(Right.a,Right.b,tmp[4]);
_Mul_FP4(Left.Tower.Gamma,tmp[1],Result.a);
_Add_FP4(Result.a,tmp[0],Result.a);
_Mul_FP4(tmp[3],tmp[4],Result.b);
_Sub_FP4(Result.b,tmp[0],tmp[2]);
_Sub_FP4(tmp[2],tmp[1],Result.b);
end;

        {********** Multiply FP4 with FP8 *****************}
procedure _Mul_FP4_FP8(const Left: Fp4Int; Right: Fp8Int; var Result: Fp8Int);
begin
Result.Tower:=Right.Tower;
_Mul_FP4(Right.a,Left,Result.a);
_Mul_FP4(Right.b,Left,Result.b);
end;

        {********** Multiply FP2 with FP8 *****************}
procedure _Mul_FP2_FP8(const Left: Fp2Int; Right: Fp8Int; var Result: Fp8Int);
begin
Result.Tower:=Right.Tower;
_Mul_FP2_FP4(Left,Right.a,Result.a);
_Mul_FP2_FP4(Left,Right.b,Result.b);
end;

        {********** Multiply FP8 with W *****************}
procedure _Mul_FP8_By_W(const value: Fp8Int; var Result: Fp8Int);  // Value should be different than Result
begin
Result.Tower:=Value.Tower;
Result.b:=Value.a;
_Mul_FP4(Value.b,Value.Tower^.Gamma,Result.a);
end;

        {********** Get Square of an FP8 *****************}
procedure _Sqr_FP8(const value: Fp8Int; var Result: Fp8Int);
var tmp: array [0..3] of Fp4Int;
begin
Result.Tower:=Value.Tower;
      {       Using Karatsuba Squering    }
_Sqr_FP4(Value.a,tmp[0]);
_Sqr_FP4(Value.b,tmp[1]);
_Add_FP4(Value.a,Value.b,tmp[3]);
_Mul_FP4(tmp[1],Value.Tower.Gamma,Result.a);
_Add_FP4(Result.a,tmp[0],Result.a);
_Sqr_FP4(tmp[3],tmp[3]);
_Sub_FP4(tmp[3],tmp[0],tmp[3]);
_Sub_FP4(tmp[3],tmp[1],Result.b);
end;

        {**********   Compare two FP8 integers *************}
function _Equals_FP8(const Left,Right:Fp8Int):boolean;
begin
result:=(_Equals_FP4(Left.a,Right.a)) and (_Equals_FP4(Left.b,Right.b));
end;

        {********** Get inverse of an FP8 *****************}
procedure _Inv_FP8(const value:Fp8Int;var Result:Fp8Int);
var t:array[0..9] of Fp4Int;
begin
Result.Tower:=Value.Tower;
_Sqr_FP4(value.a,t[0]);
_Sqr_FP4(value.b,t[1]);
_Mul_FP4(t[1],Value.Tower.Gamma,t[2]);
_Sub_FP4(t[0],t[2],t[0]);
_Inv_FP4(t[0],t[0]);
_Mul_FP4(Value.a,t[0],Result.a);
_Mul_FP4(Value.b,t[0],t[1]);
_Neg_FP4(t[1],Result.b);
end;

        {********** Get negative of an FP8 *****************}
procedure _Neg_FP8(const Value: Fp8Int; var Result: Fp8Int);
begin
Result.Tower:=Value.Tower;
_Neg_FP4(Value.a,Result.a);
_Neg_FP4(Value.b,Result.b);
end;

        {********** Conjugate of an FP8 *****************}
procedure _Conjugate_FP8(const Source: Fp8Int; var Result: Fp8Int);
begin
Result.Tower := Source.Tower;
Result.a := Source.a;
if not Source.b.IsZero then
    _Neg_FP4(Source.b, Result.b)
  else
    Result.b := Source.b;
end;


        {********** Raise an FP8 to a LInt power ************}
        // Naive Square and Multiply
procedure _Pow_FP8(const value:Fp8Int;Exponent:LInt;var Result:Fp8Int);
var i:word;
    tmp:Fp8Int;
begin
Result.Tower:= Value.Tower;
Result.a:=Value.a;
Result.b:=Value.b;
tmp:=Value;
if (Exponent>1) then
for i:=Exponent.BitLength-2 downto 0 do begin
                                        _Sqr_FP8(Result,Result);
                                        if _Is_BitSet_At(Exponent,i) then _Mul_FP8(Result,tmp,Result);
                                        end;
end;

        {****** Convert an FP8 to a Decimal String **********}
function _FP8_To_DecimalString(const value:Fp8Int):string;
begin
if Value.b.IsZero then begin
                       if Value.a.IsZero then Result:='0'
                       else Result:=Value.a.ToDecimalString;
                       end
else begin
     if Value.b.IsOne then begin
                           if Value.a.IsZero then Result:='w'
                           else Result:='w +'+Value.a.ToDecimalString;
                           end
     else begin
          if Value.a.IsZero then Result:='('+Value.b.ToDecimalString+') * w'
          else Result:='('+Value.b.ToDecimalString+') * w +'+Value.a.ToDecimalString;
          end;
     end;
end;

        {****** Convert an FP8 to a Hexadecimal String *******}
function _FP8_To_HexString(const value:Fp8Int):string;
begin
if Value.b.IsZero then begin
                       if Value.a.IsZero then Result:='0'
                       else Result:=Value.a.toHexString;
                       end
else begin
     if Value.b.IsOne then begin
                           if Value.a.IsZero then Result:='v'
                           else Result:='v +'+Value.a.toHexString;
                           end
     else begin
          if Value.a.IsZero then Result:='('+Value.b.toHexString+') * v'
          else Result:='('+Value.b.toHexString+') * v +'+Value.a.toHexString;
          end;
     end;
end;

        {****** Convert String to an FP8 (Decimal/Hex)*******}
procedure _FP8_From_Strings(const valueA,valueB,valueC,valueD:String;var Result:Fp8Int);
begin
Result.a.a.SetFormString(valueA);
Result.a.b.SetFormString(valueB);
Result.b.a.SetFormString(valueC);
Result.b.b.SetFormString(valueD);
end;

{*******************************************************************************}
      ///  Definitions of an FP8 integer operators and functions
{*******************************************************************************}


class operator Fp8Int.Add(const Left, Right: Fp8Int): Fp8Int;
begin
_Add_FP8(Left,Right,Result);
end;

class operator Fp8Int.Subtract(const Left, Right: Fp8Int): Fp8Int;
begin
_Sub_FP8(Left,Right,Result);
end;

class operator Fp8Int.Multiply(const Left, Right: Fp8Int): Fp8Int;
begin
_Mul_FP8(Left,Right,Result);
end;

class operator Fp8Int.Multiply(const Left: Fp4Int;  Right: Fp8Int): Fp8Int;
begin
_Mul_FP4_FP8(Left,Right,Result);
end;

function Fp8Int.MultiplyByGamma: Fp8Int;
begin
_Mul_FP8_By_W(Self,Result);
end;

function Fp8Int.Sqr: Fp8Int;
begin
_Sqr_FP8(Self,Result);
end;

class operator Fp8Int.Equal(const Left, Right: Fp8Int): Boolean;
begin
result:=_Equals_FP8(Left,Right);
end;

function Fp8Int.Inverse: Fp8Int;
begin
_Inv_FP8(Self,Result);
end;

function Fp8Int.IsOne: boolean;
begin
Result:=(a.IsOne)and(b.IsZero);
end;

function Fp8Int.IsZero: boolean;
begin
Result:=(a.IsZero)and(b.IsZero);
end;

class operator Fp8Int.Negative(const Value: Fp8Int): Fp8Int;
begin
_Neg_FP8(Value,Result);
end;

class operator Fp8Int.NotEqual(const Left, Right: Fp8Int): Boolean;
begin
Result:=not _Equals_FP8(Left,Right)
end;

function Fp8Int.Pow(Exponent: LInt): Fp8Int;
begin
_Pow_FP8(Self,Exponent,Result);
end;

procedure Fp8Int.SetTowerParams(TParam:PtrTowerParams24);
begin
Tower:=TParam;
New(Tower.Tower8paprms);
Tower.Tower8paprms.Sigma:=Tower.sigma;
Tower.Tower8paprms.pmod8:=Tparam.pmod8;
Tower.Tower8paprms.FrobeniusPi3xSigmaSqr:=Tparam.FrobeniusPi3xSigmaSqr;
Tower.Tower8paprms.FieldParam:=Tower.FieldParam;
a.SetTowerParams(Tower.Tower8paprms);
b.SetTowerParams(Tower.Tower8paprms);
end;

procedure Fp8Int.SetToZero;
begin
a.a.a:=0;
a.a.b:=0;
a.b.a:=0;
a.b.b:=0;
b.a.a:=0;
b.a.b:=0;
b.b.a:=0;
b.b.b:=0;
end;

procedure Fp8Int.SetFromStrings(a, b,c,d: String);
begin
_FP8_From_Strings(a,b,c,d,Self);
end;

procedure Fp8Int.SetToOne;
begin
a.a.a:=1;
a.a.b:=0;
a.b.a:=0;
a.b.b:=0;
b.a.a:=0;
b.a.b:=0;
b.b.a:=0;
b.b.b:=0;
end;

procedure Fp8Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
end;

function Fp8Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.a.BitLength mod 8) * 8;
Setlength(Result,size*8);
Move(a.a.a.Data.i8[0],Result[0],size);
Move(a.a.b.Data.i8[0],Result[size],size);
Move(a.b.a.Data.i8[0],Result[2*size],size);
Move(a.b.b.Data.i8[0],Result[3*size],size);
Move(b.a.a.Data.i8[0],Result[6*Size],size);
Move(b.a.b.Data.i8[0],Result[7*size],size);
Move(b.b.a.Data.i8[0],Result[8*size],size);
Move(b.b.b.Data.i8[0],Result[9*size],size);
end;

function Fp8Int.toDecimalString: string;
begin
Result:=_FP8_To_DecimalString(Self)
end;

function Fp8Int.toHexString: string;
begin
Result:=_FP8_To_HexString(Self)
end;

end.
