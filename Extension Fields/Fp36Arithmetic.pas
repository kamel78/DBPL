unit Fp36Arithmetic;

interface

uses System.SysUtils, VCL.dialogs, LargeIntegers, Fp12Arithmetic, FP6Arithmetic,FP2Arithmetic,FP8Arithmetic, GeneralTypes;

{ *******************************************************************************************
  Développer par Faraoun Kamel Mohamed
  Université Djilali Liabes -Sidi Bel Abbes - Algérie
  kamel_mh@yahoo.fr

  Arithmetic computation over FP36 (finit field with primecaracteristiv p ). FP36 is  the sextic
  extention of FP6 (Tower extention of order 6 for FP6). With respect to the irreducible polynomial
  Z^2-Zeta=0  (Z^2=Zeta).  Elements Are in the form a+b*Z+c*Z^2, where a,b and c are from Fp12.
  ******************************************************************************************** }

type

  FP36Int = record
    a, b, c: Fp12Int;
    Tower: PtrTowerParams24;
  public
    class operator Add(const Left, Right: FP36Int): FP36Int;
    class operator Subtract(const Left, Right: FP36Int): FP36Int;
    class operator Multiply(const Left, Right: FP36Int): FP36Int;
    class operator Multiply(const Left: Fp12Int; Right: FP36Int): FP36Int;
    class operator Negative(const Value: FP36Int): FP36Int;
    class operator Equal(const Left, Right: FP36Int): Boolean;
    class operator NotEqual(const Left, Right: FP36Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse: FP36Int;
    function Pow(Exponent: LInt;Poweringmode:FpPoweringMode =pmNormal): FP36Int;
    function IsZero: Boolean;
    function IsOne: Boolean;
    function Sqr: FP36Int;
    function MultiplyByZetta: FP36Int;
    procedure SetToRandom;
    procedure SetToOne;
    procedure SetToZero;
    procedure SetTowerParams(TParam1: PtrTowerParams24);
  end;

  procedure _Add_FP36(const Left, Right: FP36Int; var Result: FP36Int);
  procedure _Sub_FP36(Const Left, Right: FP36Int; var Result: FP36Int);
  procedure _Mul_FP36(const Left, Right: FP36Int; var Result: FP36Int);
  procedure _Mul_Fp12_FP36(const Left: Fp12Int; Right: FP36Int;var Result: FP36Int);
  procedure _Mul_FP36_Zetta(const Value: FP36Int; var Result: FP36Int);
  procedure _Sqr_FP36(const Value: FP36Int; var Result: FP36Int);
  procedure _Conjugate_FP36(const Value: FP36Int; var Result: FP36Int);
  function  _Equals_FP36(const Left, Right: FP36Int): Boolean;
  procedure _Inv_FP36(const Value: FP36Int; var Result: FP36Int);
  procedure _Neg_FP36(const Value: FP36Int; var Result: FP36Int);
  procedure _Pow_FP36(const Value: FP36Int; Exponent: LInt; var Result: FP36Int;Poweringmode :FpPoweringMode=pmNormal);
  function  _FP36_To_DecimalString(const Value: FP36Int): string;
  function  _FP36_To_HexString(const Value: FP36Int): string;
  procedure _Pow_FP36_P_i(const Value: FP36Int; power: integer;var Result: FP36Int);
  procedure _Sparse_Mul_FP36(const Left, Right: FP36Int; var Result: FP36Int;TwistMode:TTwistModel);
  procedure _Compressed_FP36(const Value: FP36Int; var Result: FP36Int);
  procedure _DeCompressed_FP36(const Value: FP36Int; var Result: FP36Int);
  procedure _Compressed_Sqr_FP36(const Value: FP36Int; var Result: FP36Int);
  procedure _Unitary_Sqr_FP36(const Value: FP36Int; var Result: FP36Int);



implementation

{ ******************************************************************************* }
                        // Procedures for FP36 Arithmetic
{ ******************************************************************************* }
      { **********   Add two FP36 integers ***************** }
procedure _Add_FP36(const Left, Right: FP36Int; var Result: FP36Int);
begin
Result.Tower := Left.Tower;
_Add_Fp12(Left.a, Right.a, Result.a);
_Add_Fp12(Left.b, Right.b, Result.b);
_Add_Fp12(Left.c, Right.c, Result.c);
end;

      { **********   Substract two FP36 integers *********** }
procedure _Sub_FP36(Const Left, Right: FP36Int; var Result: FP36Int);
begin
Result.Tower := Left.Tower;
_Sub_Fp12(Left.a, Right.a, Result.a);
_Sub_Fp12(Left.b, Right.b, Result.b);
_Sub_Fp12(Left.c, Right.c, Result.c);
end;

      { ********** Multiply two FP36 integers ***************** }
procedure _Mul_FP36(const Left, Right: FP36Int; var Result: FP36Int);
var tmp: array [0 .. 9] of Fp12Int;
begin
Result.SetTowerParams(Left.Tower);
  { Multiplication Using Karatsuba Multiplication }
_Mul_Fp12(Left.a, Right.a, tmp[0]);
_Mul_Fp12(Left.b, Right.b, tmp[1]);
_Mul_Fp12(Left.c, Right.c, tmp[2]);
_Add_Fp12(Left.b, Left.c, tmp[3]);
_Add_Fp12(Right.b, Right.c, tmp[4]);
_Add_Fp12(Left.a, Left.b, tmp[5]);
_Add_Fp12(Right.a, Right.b, tmp[6]);
_Add_Fp12(Left.a, Left.c, tmp[7]);
_Add_Fp12(Right.a, Right.c, tmp[8]);
_Mul_Fp12(tmp[3], tmp[4], Result.a);
_Sub_Fp12(Result.a, tmp[1], Result.a);
_Sub_Fp12(Result.a, tmp[2], Result.a);
_Mul_Fp12_By_W(Result.a,Left.Tower.Gamma2, tmp[4]);
  /// value should be <> than result
_Add_Fp12(tmp[0], tmp[4], Result.a);
_Mul_Fp12(tmp[5], tmp[6], Result.b);
_Sub_Fp12(Result.b, tmp[0], Result.b);
_Sub_Fp12(Result.b, tmp[1], Result.b);
_Mul_Fp12_By_W(tmp[2], Left.Tower.Gamma2,tmp[9]);
_Add_Fp12(Result.b, tmp[9], Result.b);
_Mul_Fp12(tmp[7], tmp[8], Result.c);
_Sub_Fp12(Result.c, tmp[0], Result.c);
_Sub_Fp12(Result.c, tmp[2], Result.c);
_Add_Fp12(Result.c, tmp[1], Result.c);
end;

        { ********** Multiply two FP36 integers ***************** }
procedure _Sparse_Mul_FP36(const Left, Right: FP36Int; var Result: FP36Int;TwistMode:TTwistModel);
      ///  if Twistmode is D-type then Right.c=0 and right.b.b=0
      ///  if Twistmode is M-type then Right.b=0 and right.c.b=0
var tmp: array [0 .. 9] of Fp12Int;
begin
Result.SetTowerParams(Left.Tower);
  { Sparse Multiplication Using Karatsuba Multiplication }
if TwistMode=twDType then begin
                          _Mul_Fp12(Left.a, Right.a, tmp[0]);
                          //_Sparse_Mul_Fp12(Left.b, Right.b, tmp[1]);
                          _Mul_Fp12(Left.b, Right.b, tmp[1]);
                          _Add_Fp12(Left.b, Left.c, tmp[3]);
                          _Add_Fp12(Left.a, Left.b, tmp[5]);
                          _Add_Fp12(Right.a, Right.b, tmp[6]);
                          _Add_Fp12(Left.a, Left.c, tmp[7]);
                          //_Sparse_Mul_Fp12(tmp[3], right.b, Result.a);
                          _Mul_Fp12(tmp[3], right.b, Result.a);
                          _Sub_Fp12(Result.a, tmp[1], Result.a);
                          _Mul_Fp12_By_W(Result.a,Result.Tower.Gamma2, tmp[4]);
                            /// value should be <> than result
                          _Add_Fp12(tmp[0], tmp[4], Result.a);
                          _Mul_Fp12(tmp[5], tmp[6], Result.b);
                          _Sub_Fp12(Result.b, tmp[0], Result.b);
                          _Sub_Fp12(Result.b, tmp[1], Result.b);
                          _Mul_Fp12(tmp[7], Right.a, Result.c);
                          _Sub_Fp12(Result.c, tmp[0], Result.c);
                          _Add_Fp12(Result.c, tmp[1], Result.c);
                          end
else begin
     Result.Tower := Left.Tower;
        { Multiplication Using Karatsuba Multiplication }
     _Mul_Fp12(Left.a, Right.a, tmp[0]);
     //_Sparse_Mul_Fp12(Left.c, Right.c, tmp[2]);
     _Mul_Fp12(Left.c, Right.c, tmp[2]);
     _Add_Fp12(Left.b, Left.c, tmp[3]);
     _Add_Fp12(Left.a, Left.b, tmp[5]);
     _Add_Fp12(Left.a, Left.c, tmp[7]);
     _Add_Fp12(Right.a, Right.c, tmp[8]);
     //_Sparse_Mul_Fp12(tmp[3], right.c, Result.a);
     _Mul_Fp12(tmp[3], right.c, Result.a);
     _Sub_Fp12(Result.a, tmp[2], Result.a);
     _Mul_Fp12_By_W(Result.a,Result.Tower.Gamma2, tmp[4]);
        /// value should be <> than result
     _Add_Fp12(tmp[0], tmp[4], Result.a);
     _Mul_Fp12(tmp[5], right.a, Result.b);
     _Sub_Fp12(Result.b, tmp[0], Result.b);
     _Mul_Fp12_By_W(tmp[2],Result.Tower.Gamma2, tmp[9]);
     _Add_Fp12(Result.b, tmp[9], Result.b);
     _Mul_Fp12(tmp[7], tmp[8], Result.c);
     _Sub_Fp12(Result.c, tmp[0], Result.c);
     _Sub_Fp12(Result.c, tmp[2], Result.c);
     end;
end;

        { ********** Multiply Fp12 with FP36 ***************** }
procedure _Mul_Fp12_FP36(const Left: Fp12Int; Right: FP36Int; var Result: FP36Int);
begin
Result.SetTowerParams(Right.Tower);
_Mul_Fp12(Right.a, Left, Result.a);
_Mul_Fp12(Right.b, Left, Result.b);
_Mul_Fp12(Right.c, Left, Result.c);
end;

        { ********** Multiply FP36 with Zetta ***************** }
procedure _Mul_FP36_Zetta(const Value: FP36Int; var Result: FP36Int);
// Value should be different than Result
begin
Result.SetTowerParams(Value.Tower);
Result.b := Value.a;
Result.c := Value.b;
_Mul_Fp12_By_W(Value.c,Result.Tower.Gamma2, Result.a);
end;

        { ********** Get Compressed Square of an FP36 (Karabina Approach)******* }
/// Squaring into the Cyclotmoic sub-groupe
procedure _Compressed_FP36(const Value: FP36Int; var Result: FP36Int);
var t: array [0 .. 5] of FP6Int;
begin
Result.SetTowerParams(Value.Tower);
Result := Value;
Result.a.a.SetToZero;
Result.a.b.SetToZero;
end;

        { ********** Get Compressed Square of an FP36 (Karabina Approach)***************** }
        /// Squaring into the Cyclotmoic sub-groupe
procedure _Compressed_Sqr_FP36(const Value: FP36Int; var Result: FP36Int);
var S4,S5,S2,S3,S45,S23:FP6int;
    t:array[0..1]of FP6Int;
begin
Result.SetTowerParams(Value.Tower);
_Add_FP6(Value.c.a,Value.c.b,S45);
_Sqr_FP6(S45,S45);
_Sqr_FP6(Value.c.a,S4);
_Sqr_FP6(Value.c.b,S5);
_Add_FP6(Value.b.a,Value.b.b,S23);
_Sqr_FP6(S23,S23);
_Sqr_FP6(Value.b.a,S2);
_Sqr_FP6(Value.b.b,S3);
_Add_FP6(Value.b.a,Value.b.a,t[1]);
_Sub_FP6(S45,S4,t[0]);
_Sub_FP6(t[0],S5,t[0]);
_Mul_FP6_By_V2(t[0],t[0]);
_Add_FP6(t[0],t[0],Result.b.a);
_Add_FP6(t[0],Result.b.a,Result.b.a);
_Add_FP6(t[1],Result.b.a,Result.b.a);
_Mul_FP6_By_V2(S5,t[0]);
_Add_FP6(t[0],S4,t[0]);
_Add_FP6(t[0],t[0],t[1]);
_Add_FP6(t[1],t[0],t[1]);
_Nomod_Sub_FP6(t[1],Value.b.b,t[1]);
_Sub_FP6(t[1],Value.b.b,Result.b.b);
_Mul_FP6_By_V2(S3,t[0]);
_Add_FP6(S2,t[0],t[0]);
_Add_FP6(t[0],t[0],t[1]);
_Add_FP6(t[1],t[0],t[1]);
_Sub_FP6(t[1],Value.c.a,t[1]);
_Sub_FP6(t[1],Value.c.a,Result.c.a);
_Sub_FP6(S23,S2,t[0]);
_Sub_FP6(t[0],S3,t[0]);
_Add_FP6(t[0],t[0],t[1]);
_Add_FP6(t[1],t[0],t[1]);
_Add_FP6(t[1],Value.c.b,t[1]);
_Add_FP6(t[1],Value.c.b,Result.c.b);
end;

        { ********** Get Compressed Square of an FP36 (Karabina Approach)***************** }
    /// Squaring into the Cyclotmoic sub-groupe with less modulation use (Faster)
{procedure _Nomod_Compressed_Sqr_FP36(const Value: FP36Int; var Result: FP36Int);
var S4,S5,S2,S3,S45,S23:FP6int;
    t:array[0..1]of FP6Int;
begin
Result.Tower := Value.Tower;
_Add_FP6(Value.c.a,Value.c.b,S45);
_Sqr_FP6(S45,S45);
_Nomod_Sqr_FP6(Value.c.a,S4);
_Nomod_Sqr_FP6(Value.c.b,S5);
_Nomod_Add_FP6(Value.b.a,Value.b.b,S23);
_Nomod_Sqr_FP6(S23,S23);
_Nomod_Sqr_FP6(Value.b.a,S2);
_Nomod_Sqr_FP6(Value.b.b,S3);
_Nomod_Add_FP6(Value.b.a,Value.b.a,t[1]);
_Nomod_Sub_FP6(S45,S4,t[0]);
_Nomod_Sub_FP6(t[0],S5,t[0]);
_Nomod_Mul_FP6_By_V(t[0],t[0]);
_Nomod_Add_FP6(t[0],t[0],Result.b.a);
_Nomod_Add_FP6(t[0],Result.b.a,Result.b.a);
_Add_FP6(t[1],Result.b.a,Result.b.a);
_Nomod_Mul_FP6_By_V(S5,t[0]);
_Nomod_Add_FP6(t[0],S4,t[0]);
_Nomod_Add_FP6(t[0],t[0],t[1]);
_Nomod_Add_FP6(t[1],t[0],t[1]);
_Nomod_Sub_FP6(t[1],Value.b.b,t[1]);
_Sub_FP6(t[1],Value.b.b,Result.b.b);
_Mul_FP6_By_V(S3,t[0]);
_Nomod_Add_FP6(S2,t[0],t[0]);
_Nomod_Add_FP6(t[0],t[0],t[1]);
_Nomod_Add_FP6(t[1],t[0],t[1]);
_Nomod_Sub_FP6(t[1],Value.c.a,t[1]);
_Sub_FP6(t[1],Value.c.a,Result.c.a);
_Nomod_Sub_FP6(S23,S2,t[0]);
_Nomod_Sub_FP6(t[0],S3,t[0]);
_Nomod_Add_FP6(t[0],t[0],t[1]);
_Nomod_Add_FP6(t[1],t[0],t[1]);
_Nomod_Add_FP6(t[1],Value.c.b,t[1]);
_Add_FP6(t[1],Value.c.b,Result.c.b);
end; }

        { ********** Get Decompressed Square of an FP36 (Karabina Approach)***************** }
        /// Squaring Compressed elements from the Cyclotmoic sub-groupe
procedure _DeCompressed_FP36(const Value: FP36Int; var Result: FP36Int);
var t: array [0 .. 3] of FP6Int;
begin
Result.SetTowerParams(Value.Tower);
Result := Value;
if Value.b.a.IsZero then begin
                         _Mul_FP6(Value.c.a,Value.c.b,t[0]);
                         _Add_FP6(t[0],t[0],t[0]);
                         _Inv_FP6(Value.b.b,t[1]);
                         _Mul_FP6(t[0],t[1],Result.a.b);
                         _Sqr_FP6(Result.a.b,t[0]);
                         _Add_FP6(t[0],t[0],t[0]);
                         _Mul_FP6(Value.b.b,Value.c.a,t[1]);
                         _Mul_FP_FP6(Lint(3),t[1],t[1]);
                         _Sub_FP6(t[0],t[1],t[0]);
                         _Mul_FP6(t[0],Value.Tower.Gamma2,Result.a.a);
                         _Inc_LInt(Result.a.a.a.a,1);
                         end
else  begin
      _Sqr_FP6(Value.c.b, t[0]);
      _Mul_FP6_By_V2(t[0],t[0]);
      _Sqr_FP6(Value.c.a, t[1]);
      _Add_FP6(t[1], t[1], t[2]);
      _Add_FP6(t[1], t[2], t[2]);    /// 3*t[1]
      _Add_FP6(t[0], t[2], t[2]);
      _Sub_FP6(t[2],Value.b.b,t[2]);
      _Sub_FP6(t[2],Value.b.b,t[2]);
      _Add_FP6(Value.b.a,Value.b.a,t[3]);
      _Add_FP6(t[3],t[3],t[3]);
      _Inv_FP6(t[3],t[1]);
      _Mul_FP6(t[2],t[1],Result.a.b);
      _Sqr_FP6(Result.a.b,t[2]);
      _Add_FP6(t[2],t[2],t[0]);
      _Mul_FP6(Value.b.a,Value.c.b,t[1]);
      _Add_FP6(t[0],t[1],t[0]);
      _Mul_FP6(Value.b.b,Value.c.a,t[1]);
      _Sub_FP6(t[0],t[1],t[0]);
      _Sub_FP6(t[0],t[1],t[0]);
      _Sub_FP6(t[0],t[1],t[0]);
      _Mul_FP6_By_V2(t[0],Result.a.a);
      _Inc_LInt(Result.a.a.a.a,1);
      end;
end;

        { ********** Get Square of an FP36 , Only in the Cyclotomic Group ***************** }
procedure _Unitary_Sqr_FP36(const Value: FP36Int; var Result: FP36Int);
        // Granger & Scott PKC 2010
var tmp: array [0 .. 3] of Fp12Int;
begin
Result.SetTowerParams(Value.Tower);
_Sqr_Fp12(Value.a,tmp[0]);
_Add_Fp12(tmp[0],tmp[0],tmp[1]);
_Add_Fp12(tmp[1],tmp[0],tmp[1]);
_Conjugate_Fp12(Value.a,tmp[0]);
_Add_Fp12(tmp[0],tmp[0],tmp[0]);
_Sub_Fp12(tmp[1],tmp[0],Result.a);
_Sqr_Fp12(Value.c,tmp[0]);
_Mul_Fp12_By_W(tmp[0],Result.Tower.Gamma2, tmp[1]);
_Add_Fp12(tmp[1],tmp[1],tmp[0]);
_Add_Fp12(tmp[1],tmp[0],tmp[1]);
_Sqr_Fp12(Value.b,tmp[0]);
_Add_Fp12(tmp[0],tmp[0],tmp[2]);
_Add_Fp12(tmp[2],tmp[0],tmp[2]);
_Conjugate_Fp12(Value.b,Result.b);
_Add_Fp12(Result.b,Result.b,Result.b);
_Conjugate_Fp12(Value.c,Result.c);
_Add_Fp12(Result.c,Result.c,Result.c);
_Neg_Fp12(Result.c,Result.c);
_Add_Fp12(Result.b,tmp[1],Result.b);
_Add_Fp12(Result.c,tmp[2],Result.c);
end;

        { ********** Get Square of an FP36 ***************** }
procedure _Sqr_FP36(const Value: FP36Int; var Result: FP36Int);
var tmp: array [0 .. 5] of Fp12Int;
begin
Result.SetTowerParams(Value.Tower);
  { Using Karatsuba Squering }
_Sqr_Fp12(Value.a, tmp[0]);
_Sqr_Fp12(Value.b, tmp[1]);
_Sqr_Fp12(Value.c, tmp[2]);
_Add_Fp12(Value.a, Value.b, tmp[3]);
_Add_Fp12(Value.a, Value.c, tmp[4]);
_Add_Fp12(Value.b, Value.c, Result.a);
_Sqr_Fp12(Result.a, Result.a);
_Sub_Fp12(Result.a, tmp[1], Result.a);
_Sub_Fp12(Result.a, tmp[2], Result.a);
_Mul_Fp12_By_W(Result.a,Result.Tower.Gamma2, tmp[5]);
_Add_Fp12(tmp[5], tmp[0], Result.a);
_Sqr_Fp12(tmp[3], tmp[3]);
_Sub_Fp12(tmp[3], tmp[0], tmp[3]);
_Sub_Fp12(tmp[3], tmp[1], tmp[3]);
_Mul_Fp12_By_W(tmp[2],Result.Tower.Gamma2, Result.b);
_Add_Fp12(Result.b, tmp[3], Result.b);
_Sqr_Fp12(tmp[4], Result.c);
_Sub_Fp12(Result.c, tmp[0], Result.c);
_Sub_Fp12(Result.c, tmp[2], Result.c);
_Add_Fp12(Result.c, tmp[1], Result.c);
end;

        { **********   Compare two FP36 integers ************* }
function _Equals_FP36(const Left, Right: FP36Int): Boolean;
begin
Result := (_Equals_Fp12(Left.a, Right.a)) and (_Equals_Fp12(Left.b, Right.b))and (_Equals_Fp12(Left.c, Right.c));
end;



        { ********** Raise an FP36 to " p^i " power ************ }
procedure _Pow_FP36_P_i(const Value: FP36Int; power: integer;var Result: FP36Int);
var i: integer;
begin
Result.Tower := Value.Tower;
Result.a.SetTowerParams(Value.Tower.Tower8paprms);
Result.b.SetTowerParams(Value.Tower.Tower8paprms);
Result.c.SetTowerParams(Value.Tower.Tower8paprms);
Result := Value;
for i := 0 to power - 1 do begin
                           Result.a.a:=Result.a.a.toPowerP;
                           _Mul_FP6(Result.a.b.toPowerP,Value.Tower^.FrobeniusP_Const_Fp6[2],Result.a.b);
                           _Mul_FP6(Result.c.a.toPowerP,Value.Tower^.FrobeniusP_Const_Fp6[1],Result.c.a);
                           _Mul_FP6(Result.b.a.toPowerP,Value.Tower^.FrobeniusP_Const_Fp6[0],Result.b.a);
                           _Mul_FP6(Result.b.b.toPowerP,Value.Tower^.FrobeniusP_Const_Fp6[3],Result.b.b);
                           _Mul_FP6(Result.c.b.toPowerP,Value.Tower^.FrobeniusP_Const_Fp6[4],Result.c.b);
                           end;

end;

        { ********** Get inverse of an FP36 ***************** }
procedure _Inv_FP36(const Value: FP36Int; var Result: FP36Int);
var t: array [0 .. 9] of Fp12Int;
begin
Result.Tower := Value.Tower;
_Sqr_Fp12(Value.a, t[0]);
_Sqr_Fp12(Value.b, t[1]);
_Sqr_Fp12(Value.c, t[2]);
_Mul_Fp12(Value.a, Value.b, t[3]);
_Mul_Fp12(Value.a, Value.c, t[4]);
_Mul_Fp12(Value.b, Value.c, t[5]);
_Mul_Fp12_By_W(t[5],Result.Tower.Gamma2, t[7]);
_Sub_Fp12(t[0], t[7], t[7]);
_Mul_Fp12_By_W(t[2], Result.Tower.Gamma2, t[8]);
_Sub_Fp12(t[8], t[3], t[8]);
_Sub_Fp12(t[1], t[4], t[9]);
_Mul_Fp12(Value.a, t[7], t[6]);
_Mul_Fp12_By_W(Value.c, Result.Tower.Gamma2, t[0]);
_Mul_Fp12(t[0], t[8], t[0]);
_Add_Fp12(t[6], t[0], t[6]);
_Mul_Fp12_By_W(Value.b, Result.Tower.Gamma2, t[0]);
_Mul_Fp12(t[0], t[9], t[0]);
_Add_Fp12(t[6], t[0], t[6]);
_Inv_Fp12(t[6], t[6]);
_Mul_Fp12(t[7], t[6], Result.a);
_Mul_Fp12(t[8], t[6], Result.b);
_Mul_Fp12(t[9], t[6], Result.c);
end;
        { ********** Get inverse of an FP36 ***************** }
{procedure _Inv_FP36(const Value: FP36Int; var Result: FP36Int);
var t: array [0 .. 5] of Fp12Int;
begin
Result.Tower := Value.Tower;
_Sqr_Fp12(Value.a, t[0]);
//_Sqr_Fp12(Value.b, t[1]);
_Sqr_Fp12(Value.c, t[2]);
_Mul_Fp12(Value.a, Value.b, t[3]);
_Mul_Fp12(Value.a, Value.c, t[4]);
_Mul_Fp12(Value.b, Value.c, t[5]);
_Mul_Fp12_By_W(t[5],Result.Tower.Gamma2, t[1]);
_Sub_Fp12(t[0], t[1], t[1]);
_Mul_Fp12_By_W(t[2], Result.Tower.Gamma2, t[5]);
_Sub_Fp12(t[5], t[3], t[5]);
_Sqr_Fp12(Value.b, t[2]);
_Sub_Fp12(t[2], t[4], t[4]);
_Mul_Fp12(Value.a, t[1], t[2]);
_Mul_Fp12_By_W(Value.c, Result.Tower.Gamma2, t[0]);
_Mul_Fp12(t[0], t[5], t[0]);
_Add_Fp12(t[2], t[0], t[2]);
_Mul_Fp12_By_W(Value.b, Result.Tower.Gamma2, t[0]);
_Mul_Fp12(t[0], t[4], t[0]);
_Add_Fp12(t[2], t[0], t[2]);
_Inv_Fp12(t[2], t[2]);
_Mul_Fp12(t[1], t[2], Result.a);
_Mul_Fp12(t[5], t[2], Result.b);
_Mul_Fp12(t[4], t[2], Result.c);
end;}

        { ********** Get negative of an FP36 ***************** }
procedure _Neg_FP36(const Value: FP36Int; var Result: FP36Int);
begin
Result.Tower := Value.Tower;
_Neg_Fp12(Value.a, Result.a);
_Neg_Fp12(Value.b, Result.b);
_Neg_Fp12(Value.c, Result.c);
end;

        { ********** Conjugate an FP36 ***************** }
procedure _Conjugate_FP36(const Value: FP36Int; var Result: FP36Int);
begin
  Result.Tower := Value.Tower;
  _Conjugate_Fp12(Value.a, Result.a);
  _Conjugate_Fp12(Value.b, Result.b);
  _Neg_Fp12(Result.b, Result.b);
  _Conjugate_Fp12(Value.c, Result.c);
// Faster Conjugate of FP36
{if Value<>Result then Result:=value;
_Neg_FP2(Result.a.b.a,Result.a.b.a);
_Neg_FP2(Result.a.b.b,Result.a.b.b);
_Neg_FP2(Result.b.a.a,Result.b.a.a);
_Neg_FP2(Result.b.a.b,Result.b.a.b);
_Neg_FP2(Result.c.b.a,Result.c.b.a);
_Neg_FP2(Result.c.b.b,Result.c.b.b);}
end;

        { ********** Raise an FP36 to a LInt power ************ }
procedure _Pow_FP36(const Value: FP36Int; Exponent: LInt; var Result: FP36Int;Poweringmode :FpPoweringMode=pmNormal);
var  i: word;
     tmp,zi: FP36Int;
     nafexpo,t1,t2:LIntArrayForm;
     SizeNaf,Sizebin:Integer;
begin
Result.Tower := Value.Tower;
case Poweringmode of
  pmNormal: begin
              Result:=Value;
              tmp := Value;
              if (Exponent.Absolute > 1) then
                for i := Exponent.BitLength - 2 downto 0 do begin
                                                            _Sqr_FP36(Result, Result);
                                                            if _Is_BitSet_At(Exponent, i) then _Mul_FP36(Result, tmp, Result);
                                                            end;
              end;
  pmKarbina:begin
               //  Raise an FP36 to a LInt power (Karabina Approach)************
               //  Value should verify the condition (Value)^(p^36+1)=1 :Cyclotomic sub-group
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
                                     if NafExpo[0]=-1 then _Conjugate_FP36(Result,Result);
                                     end
              else begin
                   Result.SetToOne;
                   Result.SetTowerParams(Value.Tower);
                   end;
              for i:=1 to Length(NafExpo)-1 do begin
                                               _Compressed_Sqr_FP36(tmp,tmp);
                                               if nafexpo[i]=1 then begin
                                                                    _DeCompressed_FP36(tmp,zi);
                                                                    _Mul_FP36(zi,Result,Result);
                                                                    end
                                               else if nafexpo[i]=-1 then begin
                                                                          _DeCompressed_FP36(tmp,zi);
                                                                          _Conjugate_FP36(zi,zi);
                                                                          _Mul_FP36(zi,Result,Result);
                                                                          end
                                               end;
              end;
  end;
if _IsNeg(Exponent) then _Conjugate_FP36(Result,Result);
end;

        { ****** Convert an FP36 to a Decimal String ********** }
function _FP36_To_DecimalString(const Value: FP36Int): string;
begin
if Value.c.IsZero then begin
    if Value.b.IsZero then begin
                           if Value.a.IsZero then
                           Result := '0'
                           else Result := Value.a.ToDecimalString;
                           end
    else begin
         if Value.b.IsOne then begin
            if Value.a.IsZero then Result := 'z'
            else Result := 'z +' + Value.a.ToDecimalString;
         end
         else  begin
              if Value.a.IsZero then Result := '(' + Value.b.ToDecimalString + ') * z'
              else Result := '(' + Value.b.ToDecimalString + ') * z +' +Value.a.ToDecimalString;
        end;
    end
  end
  else begin
    if Value.c.IsOne then begin
                          if Value.b.IsZero then begin
                          if Value.a.IsZero then Result := 'z^2'
                          else Result := 'z^2 + ' + Value.a.ToDecimalString;
                          end
      else begin
        if Value.b.IsOne then begin
                              if Value.a.IsZero then Result := 'z^2 + z'
                              else Result := 'z^2 + z +' + Value.a.ToDecimalString;
                              end
        else begin
             if Value.a.IsZero then Result := 'z^2 +(' + Value.b.ToDecimalString + ') * z'
             else  Result := 'z^2 + (' + Value.b.ToDecimalString + ') * z +' +Value.a.ToDecimalString;
             end;
        end
    end
    else begin
      if Value.b.IsZero then begin
                             if Value.a.IsZero then Result := '(' + Value.c.ToDecimalString + ') * z^2'
                             else Result := '(' + Value.c.ToDecimalString + ') * z^2 + ' +Value.a.ToDecimalString;
                             end
      else begin
           if Value.b.IsOne then begin
           if Value.a.IsZero then Result := '(' + Value.c.ToDecimalString + ') * z^2 +' + 'z'
           else Result := '(' + Value.c.ToDecimalString + ') * z^2 +' + ' z +' +Value.a.ToDecimalString;
        end
        else begin
             if Value.a.IsZero then Result := '(' + Value.c.ToDecimalString + ') * z^2 +' + '(' +Value.b.ToDecimalString + ') * z'
             else Result := '(' + Value.c.ToDecimalString + ') * z^2 +' + '(' +Value.b.ToDecimalString + ') * z +' + Value.a.ToDecimalString;
             end;
          end
      end;
end;
end;

        { ****** Convert an FP36 to a Hexadecimal String ******* }
function _FP36_To_HexString(const Value: FP36Int): string;
begin
  if Value.c.IsZero then begin
                         if Value.b.IsZero then begin
                                                if Value.a.IsZero then Result := '0'
                                                else Result := Value.a.ToHexString;
                                                end
                         else begin
                              if Value.b.IsOne then begin
                                                    if Value.a.IsZero then Result := 'z'
                                                    else Result := 'z +' + Value.a.ToHexString;
                                                    end
                              else begin
                                   if Value.a.IsZero then Result := '(' + Value.b.ToHexString + ') * z'
                                   else Result := '(' + Value.b.ToHexString + ') * z +' + Value.a.ToHexString;
                                   end;
                              end
                         end
  else  begin
        if Value.c.IsOne then begin
                              if Value.b.IsZero then begin
                                                     if Value.a.IsZero then Result := 'z^2'
                                                     else Result := 'z^2 + ' + Value.a.ToHexString;
                                                     end
                              else begin
                                   if Value.b.IsOne then begin
                                                         if Value.a.IsZero then Result := 'z^2 + z'
                                                         else Result := 'z^2 + z +' + Value.a.ToHexString;
                                                         end
                                   else begin
                                        if Value.a.IsZero then Result := 'z^2 +(' + Value.b.ToHexString + ') * z'
                                        else Result := 'z^2 + (' + Value.b.ToHexString + ') * z +' +Value.a.ToHexString;
                                        end;
                                   end
                              end
        else begin
             if Value.b.IsZero then begin
                                    if Value.a.IsZero then Result := '(' + Value.c.ToHexString + ') * z^2'
                                    else Result := '(' + Value.c.ToHexString + ') * z^2 + ' +Value.a.ToHexString;
                                    end
             else begin
                  if Value.b.IsOne then begin
                                        if Value.a.IsZero then Result := '(' + Value.c.ToHexString + ') * z^2 +' + 'z'
                                        else Result := '(' + Value.c.ToHexString + ') * z^2 +' + ' z +' +
                                        Value.a.ToHexString;
                                        end
                  else begin
                       if Value.a.IsZero then Result := '(' + Value.c.ToHexString + ') * z^2 +' + '(' +Value.b.ToHexString + ') * z'
                       else Result := '(' + Value.c.ToHexString + ') * z^2 +' + '(' +Value.b.ToHexString + ') * z +' + Value.a.ToHexString;
                       end;
                  end
             end;
        end;
end;

{ ******************************************************************************* }
/// Definitions of an FP36 integer operators and functions
{ ******************************************************************************* }

class operator FP36Int.Add(const Left, Right: FP36Int): FP36Int;
begin
  _Add_FP36(Left, Right, Result);
end;

class operator FP36Int.Subtract(const Left, Right: FP36Int): FP36Int;
begin
  _Sub_FP36(Left, Right, Result);
end;

class operator FP36Int.Multiply(const Left, Right: FP36Int): FP36Int;
begin
  _Mul_FP36(Left, Right, Result);
end;

class operator FP36Int.Multiply(const Left: Fp12Int; Right: FP36Int): FP36Int;
begin
  _Mul_Fp12_FP36(Left, Right, Result);
end;

function FP36Int.MultiplyByZetta: FP36Int;
begin
  _Mul_FP36_Zetta(Self, Result);
end;

function FP36Int.Sqr: FP36Int;
begin
  _Sqr_FP36(Self, Result);
end;

class operator FP36Int.Equal(const Left, Right: FP36Int): Boolean;
begin
  Result := _Equals_FP36(Left, Right);
end;

function FP36Int.Inverse: FP36Int;
begin
  _Inv_FP36(Self, Result);
end;

function FP36Int.IsOne: Boolean;
begin
  Result := (a.IsOne) and (b.IsZero) and (c.IsZero);
end;

function FP36Int.IsZero: Boolean;
begin
  Result := (a.IsZero) and (b.IsZero) and (c.IsZero);
end;

class operator FP36Int.Negative(const Value: FP36Int): FP36Int;
begin
  _Neg_FP36(Value, Result);
end;

class operator FP36Int.NotEqual(const Left, Right: FP36Int): Boolean;
begin
  Result := not _Equals_FP36(Left, Right)
end;

function FP36Int.Pow(Exponent: LInt;Poweringmode:FpPoweringMode=pmNormal): FP36Int;
begin
  _Pow_FP36(Self, Exponent, Result,Poweringmode);
end;

procedure FP36Int.SetTowerParams(TParam1: PtrTowerParams24);
begin
  Tower := TParam1;
  a.SetTowerParams(TParam1.Tower8paprms);
  b.SetTowerParams(TParam1.Tower8paprms);
  c.SetTowerParams(TParam1.Tower8paprms);
end;

procedure FP36Int.SetToZero;
begin
  a.SetToOne;
  b.SetToZero;
  c.SetToZero;
end;

procedure FP36Int.SetToOne;
begin
  a.SetToOne;
  b.SetToZero;
  c.SetToZero;
end;

procedure FP36Int.SetToRandom;
begin
  a.SetToRandom;
  b.SetToRandom;
  c.SetToRandom;
end;

function FP36Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.a.a.BitLength mod 8) * 8;
Setlength(Result,size*36);
Move(a.a.a.a.Data.i8[0],Result[0],size);
Move(a.a.a.b.Data.i8[0],Result[size],size);
Move(a.a.b.a.Data.i8[0],Result[2*size],size);
Move(a.a.b.b.Data.i8[0],Result[3*size],size);
Move(a.a.c.a.Data.i8[0],Result[4*size],size);
Move(a.a.c.b.Data.i8[0],Result[5*size],size);
Move(a.b.a.a.Data.i8[0],Result[6*size],size);
Move(a.b.a.b.Data.i8[0],Result[7*size],size);
Move(a.b.b.a.Data.i8[0],Result[8*size],size);
Move(a.b.b.b.Data.i8[0],Result[9*size],size);
Move(a.b.c.a.Data.i8[0],Result[10*size],size);
Move(a.b.c.b.Data.i8[0],Result[11*size],size);

Move(b.a.a.a.Data.i8[0],Result[12*size],size);
Move(b.a.a.b.Data.i8[0],Result[13*size],size);
Move(b.a.b.a.Data.i8[0],Result[14*size],size);
Move(b.a.b.b.Data.i8[0],Result[15*size],size);
Move(b.a.c.a.Data.i8[0],Result[16*size],size);
Move(b.a.c.b.Data.i8[0],Result[17*size],size);
Move(b.b.a.a.Data.i8[0],Result[18*size],size);
Move(b.b.a.b.Data.i8[0],Result[19*size],size);
Move(b.b.b.a.Data.i8[0],Result[20*size],size);
Move(b.b.b.b.Data.i8[0],Result[21*size],size);
Move(b.b.c.a.Data.i8[0],Result[22*size],size);
Move(b.b.c.b.Data.i8[0],Result[23*size],size);

Move(c.a.a.a.Data.i8[0],Result[24*size],size);
Move(c.a.a.b.Data.i8[0],Result[25*size],size);
Move(c.a.b.a.Data.i8[0],Result[26*size],size);
Move(c.a.b.b.Data.i8[0],Result[27*size],size);
Move(c.a.c.a.Data.i8[0],Result[28*size],size);
Move(c.a.c.b.Data.i8[0],Result[29*size],size);
Move(c.b.a.a.Data.i8[0],Result[30*size],size);
Move(c.b.a.b.Data.i8[0],Result[31*size],size);
Move(c.b.b.a.Data.i8[0],Result[32*size],size);
Move(c.b.b.b.Data.i8[0],Result[33*size],size);
Move(c.b.c.a.Data.i8[0],Result[34*size],size);
Move(c.b.c.b.Data.i8[0],Result[35*size],size);
end;

function FP36Int.ToDecimalString: string;
begin
  Result := _FP36_To_DecimalString(Self)
end;

function FP36Int.ToHexString: string;
begin
  Result := _FP36_To_HexString(Self)
end;

end.
