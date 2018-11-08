unit FP4Arithmetic;

interface

uses System.SysUtils, VCL.dialogs, LargeIntegers, Fp2Arithmetic, Fp6Arithmetic;

{ *******************************************************************************************
  Développer par Faraoun Kamel Mohamed
  Université Djilali Liabes -Sidi Bel Abbes - Algérie
  kamel_mh@yahoo.fr

  Arithmetic computation over Fp4 . Fp4 is  the quadratic extention of Fp2 (Tower extention of order 2 for FP2)
  with respect to the irreducible polynômial v²-Sigma (v^2=Sigma). Elements are on the form a+v.b (a,b from Fp2).

  ******************************************************************************************** }

type

  Fp4Int = record
    a, b: Fp2Int;
    Tower: PtrTowerParams12;    // Parametres of the Tower Extension scheme
  public
    class operator Add(const Left, Right: Fp4Int): Fp4Int;
    class operator Subtract(const Left, Right: Fp4Int): Fp4Int;
    class operator Multiply(const Left, Right: Fp4Int): Fp4Int;
    class operator Multiply(const Left: Fp2Int; Right: Fp4Int): Fp4Int;
    class operator Negative(const Value: Fp4Int): Fp4Int;
    class operator Equal(const Left, Right: Fp4Int): Boolean;
    class operator NotEqual(const Left, Right: Fp4Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse: Fp4Int;
    function Pow(Exponent: LInt): Fp4Int;
    function toPowerP:Fp4Int;
    function Conjugate: Fp4Int;
    function IsZero: Boolean;
    function IsOne: Boolean;
    function IsASquare: Boolean;
    function Sqr: Fp4Int;
    function MultiplyByGamma: Fp4Int;
    procedure SetToRandom;
    procedure SetToZero;
    procedure SetTowerParams(TParam: PtrTowerParams12);
    procedure SetFromStrings(a, b: AnsiString);
    procedure SetFromString(a: AnsiString);
  end;

      procedure _Add_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
      procedure _Sub_FP4(Const Left, Right: Fp4Int; var Result: Fp4Int);
      procedure _Mul_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
      procedure _Mul_FP2_FP4(const Left: Fp2Int; Right: Fp4Int; var Result: Fp4Int);
      procedure _Mul_FP_FP4(const Left: LInt; Right: Fp4Int; var Result: Fp4Int);
      procedure _Mul_FP4_By_V(const Value: Fp4Int; var Result: Fp4Int);
      function _Is_FP4_ASquare(const Value: Fp4Int): Boolean;
      procedure _Sqrt_FP4(const Value: Fp4Int; var Result: Fp4Int);
      procedure _Sqr_FP4(const Value: Fp4Int; var Result: Fp4Int);
      function _Equals_FP4(const Left, Right: Fp4Int): Boolean;
      procedure _Inv_FP4(const Value: Fp4Int; var Result: Fp4Int);
      procedure _Neg_FP4(const Value: Fp4Int; var Result: Fp4Int);
      procedure _Conjugate_FP4(const Source: Fp4Int; var Result: Fp4Int);
      procedure _Pow_FP4(const Value: Fp4Int; Exponent: LInt; var Result: Fp4Int);
      function _FP4_To_DecimalString(const Value: Fp4Int): string;
      function _FP4_To_HexString(const Value: Fp4Int): string;
      procedure _FP4_From_Strings(const valueA, valueB: String; var Result: Fp4Int);
      procedure _Pow_FP4_P(const Value: Fp4Int; var Result: Fp4Int);
      procedure _Nomod_Add_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
      procedure _Nomod_Sub_FP4(Const Left, Right: Fp4Int; var Result: Fp4Int);
      procedure _Nomod_Sqr_FP4(const Value: Fp4Int; var Result: Fp4Int);
      procedure _Nomod_Mul_FP4_By_V(const Value: Fp4Int; var Result: Fp4Int);
      procedure _Nomod_Mul_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);

implementation


{ ******************************************************************************* }
                        // Procedures for FP4 Arithmetic
{ ******************************************************************************* }

    { **********   Add two FP4 integers ***************** }
procedure _Add_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Left.Tower;
_Add_FP2(Left.a, Right.a, Result.a);
_Add_FP2(Left.b, Right.b, Result.b);
end;

    { **********   Add two FP4 integers ***************** }
procedure _Nomod_Add_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Left.Tower;
_Nomod_Add_FP2(Left.a, Right.a, Result.a);
_Nomod_Add_FP2(Left.b, Right.b, Result.b);
end;

    { **********   Substract two FP4 integers *********** }
     // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Sub_FP4(Const Left, Right: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Left.Tower;
_Nomod_Sub_FP2(Left.a, Right.a, Result.a);
_Nomod_Sub_FP2(Left.b, Right.b, Result.b);
end;

    { **********   Substract two FP4 integers ***********}
procedure _Sub_FP4(Const Left, Right: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Left.Tower;
_Sub_FP2(Left.a, Right.a, Result.a);
_Sub_FP2(Left.b, Right.b, Result.b);
end;

    { ********** Multiply two FP4 integers ***************** }
     // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Mul_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
var tmp: array [0 .. 4] of Fp2Int;
begin
  { Using Karatsuba Multiplication }
Result.Tower := Left.Tower;
_Nomod_Mul_FP2(Left.a, Right.a, tmp[0]);
_Nomod_Mul_FP2(Left.b, Right.b, tmp[1]);
_Nomod_Add_FP2(Left.a, Left.b, tmp[3]);
_Nomod_Add_FP2(Right.a, Right.b, tmp[4]);
_Nomod_Mul_FP2(tmp[1],Left.Tower.Sigma , Result.a);
_Add_FP2(Result.a, tmp[0], Result.a);
_Mul_FP2(tmp[3], tmp[4], Result.b);
_Nomod_Sub_FP2(Result.b, tmp[0], tmp[2]);
_Sub_FP2(tmp[2], tmp[1], Result.b);
end;

    { ********** Multiply two FP4 integers ***************** }
procedure _Mul_FP4(const Left, Right: Fp4Int; var Result: Fp4Int);
var tmp: array [0 .. 4] of Fp2Int;
begin
  { Using Karatsuba Multiplication }
Result.Tower := Left.Tower;
_Mul_FP2(Left.a, Right.a, tmp[0]);
_Mul_FP2(Left.b, Right.b, tmp[1]);
_Add_FP2(Left.a, Left.b, tmp[3]);
_Add_FP2(Right.a, Right.b, tmp[4]);
_Mul_FP2(tmp[1],Left.Tower.Sigma , Result.a);
_Add_FP2(Result.a, tmp[0], Result.a);
_Mul_FP2(tmp[3], tmp[4], Result.b);
_Sub_FP2(Result.b, tmp[0], tmp[2]);
_Sub_FP2(tmp[2], tmp[1], Result.b);
end;

    { ********** Multiply FP with FP4 ***************** }
procedure _Mul_FP_FP4(const Left: LInt; Right: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Right.Tower;
_Mul_FP_FP2(Left, Right.a, Result.a);
_Mul_FP_FP2(Left, Right.b, Result.b);
end;

    { ********** Multiply FP2 with FP4 ***************** }
procedure _Mul_FP2_FP4(const Left: Fp2Int; Right: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Right.Tower;
_Mul_FP2(Right.a, Left, Result.a);
_Mul_FP2(Right.b, Left, Result.b);
end;

    { ********** Multiply FP4 with V ***************** }
     // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Mul_FP4_By_V(const Value: Fp4Int; var Result: Fp4Int);
var tmp:Fp2Int;
begin
Result.Tower := Value.Tower;
_Nomod_Mul_FP2(Value.b, Value.Tower^.Sigma, tmp);
Result.b := Value.a;
Result.a:=tmp;
end;

    { ********** Multiply FP4 with V ***************** }
procedure _Mul_FP4_By_V(const Value: Fp4Int; var Result: Fp4Int);
var tmp:Fp2Int;
begin
Result.Tower := Value.Tower;
_Mul_FP2(Value.b, Value.Tower^.Sigma, tmp);
Result.b := Value.a;
Result.a:=tmp;
end;

    { ********** Get Square of an FP4 ***************** }
     // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Sqr_FP4(const Value: Fp4Int; var Result: Fp4Int);
var  tmp: array [0 .. 2] of Fp2Int;
begin
Result.Tower := Value.Tower;
  { Using Karatsuba Squering }
_Nomod_Sqr_FP2(Value.a, tmp[0]);
_Nomod_Sqr_FP2(Value.b, tmp[1]);
_Nomod_Add_FP2(Value.a, Value.b, tmp[2]);
_Nomod_Mul_FP2(tmp[1],Value.Tower.Sigma , Result.a);
_Nomod_Add_FP2(Result.a, tmp[0], Result.a);
_Nomod_Sqr_FP2(tmp[2], tmp[2]);
_Nomod_Sub_FP2(tmp[2], tmp[0], tmp[2]);
_Nomod_Sub_FP2(tmp[2], tmp[1], Result.b);
end;

    { ********** Get Square of an FP4 ***************** }
procedure _Sqr_FP4(const Value: Fp4Int; var Result: Fp4Int);
var tmp: array [0 .. 2] of Fp2Int;
begin
Result.Tower := Value.Tower;
  { Using Karatsuba Squering }
_Sqr_FP2(Value.a, tmp[0]);
_Sqr_FP2(Value.b, tmp[1]);
_nomod_Add_FP2(Value.a, Value.b, tmp[2]);
_Mul_FP2(tmp[1],Value.Tower.Sigma , Result.a);
_Add_FP2(Result.a, tmp[0], Result.a);
_Sqr_FP2(tmp[2], tmp[2]);
_Nomod_Sub_FP2(tmp[2], tmp[0], tmp[2]);
_Sub_FP2(tmp[2], tmp[1], Result.b);
end;

    {**********   Compare two FP4 integers ************* }
function _Equals_FP4(const Left, Right: Fp4Int): Boolean;
begin
Result := (_Compare_FP2(Left.a, Right.a) = 0) and (_Compare_FP2(Left.b, Right.b) = 0);
end;

    {********** Get inverse of an FP4 ***************** }
procedure _Inv_FP4(const Value: Fp4Int; var Result: Fp4Int);
var t: array [0 .. 2] of Fp2Int;
begin
Result.Tower := Value.Tower;
_Sqr_FP2(Value.a, t[0]);
_Sqr_FP2(Value.b, t[1]);
_Mul_FP2(t[1],Value.Tower.Sigma, t[2]);
_Sub_FP2(t[0], t[2], t[0]);
_Inv_FP2(t[0], t[0]);
_Mul_FP2(Value.a, t[0], Result.a);
_Mul_FP2(Value.b, t[0], t[1]);
_Neg_FP2(t[1], Result.b);
end;

    { ****** Test if an FP4 is a Sequare ******* }
function _Is_FP4_ASquare(const Value: Fp4Int): Boolean;
var t: array [0 .. 1] of Fp2Int;
begin
_Sqr_FP2(Value.b, t[0]);
_Mul_FP2(Value.Tower.Sigma, t[0], t[1]);
_Sqr_FP2(Value.a, t[0]);
_Sub_FP2(t[0], t[1], t[0]);
Result := t[0].IsASquare;
end;

    { ********** Get Square root of an FP4 ************ }
procedure _Sqrt_FP4(const Value: Fp4Int; var Result: Fp4Int);
var t: array [0 .. 2] of Fp2Int;
    delta: Fp2Int;
begin
Result.Tower := Value.Tower;
_Sqr_FP2(Value.b, t[0]);
_Mul_FP2(t[0], Value.Tower.Sigma, t[1]);
_Sqr_FP2(Value.a, delta);
_Sub_FP2(delta, t[1], delta);
if not delta.IsASquare then raise ERangeError.Create('This FP4 element is not a Square .....')
else begin
     _Sqrt_FP2(delta,t[0]);
     _Add_FP2(Value.a, t[0], t[2]);
     _Mul_FP2(t[2], Value.Tower.inv_2, t[1]);
     if not t[1].IsASquare then begin
                                _Sub_FP2(t[1], t[0], t[1]);
                                if not t[1].IsASquare then begin
                                                           _Neg_FP2(t[1], t[1]);
                                                           if not t[1].IsASquare then _Sub_FP2(t[1], t[0], t[1]);
                                                           end;
                                end;
     Result.a:=t[1].Sqrt;
     if Result.a.IsZero then Result.b.SetFormString('0+u*0')
     else begin
          t[0]:=Result.a;
          t[0] := 2 * t[0];
          _Inv_FP2(t[0], t[1]);
          _Mul_FP2(Value.b, t[1], Result.b);
          end;
     end;
end;

    { **********   Conjugate an FP4 integer ************* }
procedure _Conjugate_FP4(const Source: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Source.Tower;
Result.a := Source.a;
if not Source.b.IsZero then _Neg_FP2(Source.b, Result.b)
else Result.b := Source.b;
end;

    { ********** Get negative of an FP4 ***************** }
procedure _Neg_FP4(const Value: Fp4Int; var Result: Fp4Int);
begin
Result.Tower := Value.Tower;
_Neg_FP2(Value.a, Result.a);
_Neg_FP2(Value.b, Result.b);
end;

    { ********** Fast Raise an FP4 to " p " power ************ }
    /// for BLS 24 frobenius computation
procedure _Pow_FP4_P(const Value: Fp4Int;  var Result: Fp4Int);
begin
Result.Tower:=Value.Tower;
_Conjugate_FP2(Value.a,Result.a);
_Conjugate_FP2(Value.b,Result.b);
_Mul_FP2(Value.Tower.FrobeniusPi3xSigmaSqr,Result.b,Result.b);
end;

    { ********** Raise an FP4 to a LInt power ************ }
procedure _Pow_FP4(const Value: Fp4Int; Exponent: LInt; var Result: Fp4Int);
var i: word;
    tmp: Fp4Int;
begin
Result.Tower := Value.Tower;
Result.a := Value.a;
Result.b := Value.b;
tmp := Value;
if (Exponent > 1) then
    for i := Exponent.BitLength - 2 downto 0 do begin
                                                _Sqr_FP4(Result, Result);
                                                if _Is_BitSet_At(Exponent,i) then _Mul_FP4(Result, tmp, Result);
                                                end;
end;

    { ****** Convert an FP4 to a Decimal String ********** }
function _FP4_To_DecimalString(const Value: Fp4Int): string;
begin
if Value.b.IsZero then begin
                       if Value.a.IsZero then Result := '0'
                       else Result := Value.a.ToDecimalString;
                       end
else  begin
      if Value.b.IsOne then begin
                            if Value.a.IsZero then
                            Result := 'v'
                            else Result := 'v +' + Value.a.ToDecimalString;
                            end
      else begin
           if Value.a.IsZero then Result := '(' + Value.b.ToDecimalString + ') * v'
           else Result := '(' + Value.b.ToDecimalString + ') * v +' +
           Value.a.ToDecimalString;
           end;
      end;
end;

    { ****** Convert an FP4 to a Hexadecimal String ******* }
function _FP4_To_HexString(const Value: Fp4Int): string;
begin
if Value.b.IsZero then begin
                       if Value.a.IsZero then Result := '0'
                       else Result := Value.a.ToHexString;
                       end
else begin
     if Value.b.IsOne then begin
                           if Value.a.IsZero then Result := 'v'
                           else Result := 'v +' + Value.a.ToHexString;
                           end
     else begin
          if Value.a.IsZero then Result := '(' + Value.b.ToHexString + ') * v'
          else Result := '(' + Value.b.ToHexString + ') * v +' + Value.a.ToHexString;
          end;
     end;
end;

    { ****** Convert String to an FP4 (Decimal/Hex)******* }
procedure _FP4_From_Strings(const valueA, valueB: String; var Result: Fp4Int);
begin
Result.a.SetFormString(valueA);
Result.b.SetFormString(valueB);
end;

    { ****** Convert String to an FP4 (Decimal/Hex)******* }
procedure _FP4_From_String(const value: String; var Result: Fp4Int);
begin
if Pos('+v*',value)=-1 then Result.b.SetFormString('0+u*0')
else begin
     Result.b.SetFormString(copy(value,pos('+v*',value)+2,length(value)));
     Result.a.SetFormString(copy(value,2,pos('+v*',value)-2));
     end;
end;


{ ******************************************************************************* }
/// Definitions of an FP4 integer operators and functions
{ ******************************************************************************* }

class operator Fp4Int.Add(const Left, Right: Fp4Int): Fp4Int;
begin
_Add_FP4(Left, Right, Result);
end;

class operator Fp4Int.Subtract(const Left, Right: Fp4Int): Fp4Int;
begin
_Sub_FP4(Left, Right, Result);
end;

class operator Fp4Int.Multiply(const Left, Right: Fp4Int): Fp4Int;
begin
_Mul_FP4(Left, Right, Result);
end;

class operator Fp4Int.Multiply(const Left: Fp2Int; Right: Fp4Int): Fp4Int;
begin
_Mul_FP2_FP4(Left, Right, Result);
end;

function Fp4Int.MultiplyByGamma: Fp4Int;
begin
_Mul_FP4_By_V(Self, Result);
end;

function Fp4Int.Sqr: Fp4Int;
begin
_Sqr_FP4(Self, Result);
end;

function Fp4Int.Conjugate: Fp4Int;
begin
_Conjugate_FP4(Self, Result);
end;

class operator Fp4Int.Equal(const Left, Right: Fp4Int): Boolean;
begin
Result := _Equals_FP4(Left, Right);
end;

function Fp4Int.Inverse: Fp4Int;
begin
_Inv_FP4(Self, Result);
end;

function Fp4Int.IsASquare: Boolean;
begin
Result := _Is_FP4_ASquare(Self);
end;

function Fp4Int.IsOne: Boolean;
begin
Result := (a.IsOne) and (b.IsZero);
end;

function Fp4Int.IsZero: Boolean;
begin
Result := (a.IsZero) and (b.IsZero);
end;

class operator Fp4Int.Negative(const Value: Fp4Int): Fp4Int;
begin
_Neg_FP4(Value, Result);
end;

class operator Fp4Int.NotEqual(const Left, Right: Fp4Int): Boolean;
begin
Result := not _Equals_FP4(Left, Right)
end;

function Fp4Int.Pow(Exponent: LInt): Fp4Int;
begin
_Pow_FP4(Self, Exponent, Result);
end;

procedure Fp4Int.SetTowerParams(TParam: PtrTowerParams12);
begin
Tower := TParam;
a.SetFieldParams(TParam.FieldParam);
b.SetFieldParams(TParam.FieldParam);
end;

procedure Fp4Int.SetToZero;
begin
a.a.Data.i16[-2]:=0;
a.b.Data.i16[-2]:=0;
b.a.Data.i16[-2]:=0;
b.b.Data.i16[-2]:=0;
end;

procedure Fp4Int.SetFromStrings(a, b: AnsiString);
begin
_FP4_From_Strings(a, b, Self);
end;

procedure Fp4Int.SetFromString(a: AnsiString);
begin
_FP4_From_String(a, Self);
end;

procedure Fp4Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
end;

function Fp4Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.BitLength mod 8) * 8;
Setlength(Result,size*4);
Move(a.a.Data.i8[0],Result[0],size);
Move(a.b.Data.i8[0],Result[size],size);
Move(b.a.Data.i8[0],Result[2*size],size);
Move(b.b.Data.i8[0],Result[3*size],size);
end;

function Fp4Int.ToDecimalString: string;
begin
Result := _FP4_To_DecimalString(Self)
end;

function Fp4Int.ToHexString: string;
begin
Result := _FP4_To_HexString(Self)
end;

function Fp4Int.toPowerP: Fp4Int;
begin
_Pow_FP4_P(Self,Result);
end;

end.
