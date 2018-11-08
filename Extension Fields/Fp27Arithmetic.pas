unit Fp27Arithmetic;

interface

uses System.SysUtils, VCL.dialogs, LargeIntegers, Fp9Arithmetic, FP3Arithmetic,GeneralTypes;

{ *******************************************************************************************
  Développer par Faraoun Kamel Mohamed
  Université Djilali Liabes -Sidi Bel Abbes - Algérie
  kamel_mh@yahoo.fr

  Arithmetic computation over Fp27 (finit field with primecaracteristiv p ). Fp27 is  the cubic
  extention of Fp3 (Tower extention of order 3 for Fp3). With respect to the irreducible polynomial
  Z^3-Zeta=0  (Z^2=Zeta).  Elements Are in the form a+b*Z+c*Z^2, where a,b and c are from Fp9.
  ******************************************************************************************** }

type


  Fp27Int = record
    a, b, c: Fp9Int;
    Tower: PtrTowerParams18;
  public
    class operator Add(const Left, Right: Fp27Int): Fp27Int;
    class operator Subtract(const Left, Right: Fp27Int): Fp27Int;
    class operator Multiply(const Left, Right: Fp27Int): Fp27Int;
    class operator Multiply(const Left: Fp9Int; Right: Fp27Int): Fp27Int;
    class operator Negative(const Value: Fp27Int): Fp27Int;
    class operator Equal(const Left, Right: Fp27Int): Boolean;
    class operator NotEqual(const Left, Right: Fp27Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse: Fp27Int;
    function Pow(Exponent: LInt;Poweringmode:FpPoweringMode =pmNormal): Fp27Int;
    function IsZero: Boolean;
    function IsOne: Boolean;
    function Sqr: Fp27Int;
    function MultiplyByZetta: Fp27Int;
    procedure SetToRandom;
    procedure SetToOne;
    procedure SetToZero;
    procedure SetTowerParams(TParam: PtrTowerParams18);
    procedure SetFromStrings(a, b, c, d, e, f, g, h, i, j, k, l: String);
  end;

  procedure _Add_Fp27(const Left, Right: Fp27Int; var Result: Fp27Int);
  procedure _Sub_Fp27(Const Left, Right: Fp27Int; var Result: Fp27Int);
  procedure _Mul_Fp27(const Left, Right: Fp27Int; var Result: Fp27Int);
  procedure _Mul_Fp9_Fp27(const Left: Fp9Int; Right: Fp27Int;var Result: Fp27Int);
  procedure _Mul_Fp27_Zetta(const Value: Fp27Int; var Result: Fp27Int);
  procedure _Sqr_Fp27(const Value: Fp27Int; var Result: Fp27Int);
//  procedure _Conjugate_Fp27(const Value: Fp27Int; var Result: Fp27Int);
  function  _Equals_Fp27(const Left, Right: Fp27Int): Boolean;
  procedure _Inv_Fp27(const Value: Fp27Int; var Result: Fp27Int);
  procedure _Neg_Fp27(const Value: Fp27Int; var Result: Fp27Int);
  procedure _Pow_Fp27(const Value: Fp27Int; Exponent: LInt; var Result: Fp27Int;Poweringmode :FpPoweringMode=pmNormal);
  function  _Fp27_To_DecimalString(const Value: Fp27Int): string;
  function  _Fp27_To_HexString(const Value: Fp27Int): string;
  procedure _Fp27_From_Strings(const valueA, valueB, valueC, valueD, valueE,ValueF, ValueG, ValueH, ValueI: String;var Result: Fp27Int);
  procedure _Pow_Fp27_P(const Value: Fp27Int; var Result: Fp27Int);
  procedure _Pow_Fp27_P9(const Value: Fp27Int; var Result: Fp27Int);
  procedure _Pow_Fp27_P_i(const Value: Fp27Int;power: integer; var Result: Fp27Int);



implementation

{ ******************************************************************************* }
                        // Procedures for Fp27 Arithmetic
{ ******************************************************************************* }
      { **********   Add two Fp27 integers ***************** }
procedure _Add_Fp27(const Left, Right: Fp27Int; var Result: Fp27Int);
begin
Result.Tower := Left.Tower;
_Add_Fp9(Left.a, Right.a, Result.a);
_Add_Fp9(Left.b, Right.b, Result.b);
_Add_Fp9(Left.c, Right.c, Result.c);
end;

      { **********   Substract two Fp27 integers *********** }
procedure _Sub_Fp27(Const Left, Right: Fp27Int; var Result: Fp27Int);
begin
Result.Tower := Left.Tower;
_Sub_Fp9(Left.a, Right.a, Result.a);
_Sub_Fp9(Left.b, Right.b, Result.b);
_Sub_Fp9(Left.c, Right.c, Result.c);
end;

      { ********** Multiply two Fp27 integers ***************** }
procedure _Mul_Fp27(const Left, Right: Fp27Int; var Result: Fp27Int);
var tmp: array [0 .. 9] of Fp9Int;
begin
Result.Tower := Left.Tower;
  { Multiplication Using Karatsuba Multiplication }
_Mul_Fp9(Left.a, Right.a, tmp[0]);
_Mul_Fp9(Left.b, Right.b, tmp[1]);
_Mul_Fp9(Left.c, Right.c, tmp[2]);
_Add_Fp9(Left.b, Left.c, tmp[3]);
_Add_Fp9(Right.b, Right.c, tmp[4]);
_Add_Fp9(Left.a, Left.b, tmp[5]);
_Add_Fp9(Right.a, Right.b, tmp[6]);
_Add_Fp9(Left.a, Left.c, tmp[7]);
_Add_Fp9(Right.a, Right.c, tmp[8]);
_Mul_Fp9(tmp[3], tmp[4], Result.a);
_Sub_Fp9(Result.a, tmp[1], Result.a);
_Sub_Fp9(Result.a, tmp[2], Result.a);
_Mul_Fp9_By_V(Result.a, tmp[4]);
  /// value should be <> than result
_Add_Fp9(tmp[0], tmp[4], Result.a);
_Mul_Fp9(tmp[5], tmp[6], Result.b);
_Sub_Fp9(Result.b, tmp[0], Result.b);
_Sub_Fp9(Result.b, tmp[1], Result.b);
_Mul_Fp9_By_V(tmp[2], tmp[9]);
_Add_Fp9(Result.b, tmp[9], Result.b);
_Mul_Fp9(tmp[7], tmp[8], Result.c);
_Sub_Fp9(Result.c, tmp[0], Result.c);
_Sub_Fp9(Result.c, tmp[2], Result.c);
_Add_Fp9(Result.c, tmp[1], Result.c);
end;


        { ********** Multiply Fp9 with Fp27 ***************** }
procedure _Mul_Fp9_Fp27(const Left: Fp9Int; Right: Fp27Int; var Result: Fp27Int);
begin
Result.Tower := Right.Tower;
_Mul_Fp9(Right.a, Left, Result.a);
_Mul_Fp9(Right.b, Left, Result.b);
_Mul_Fp9(Right.c, Left, Result.c);
end;

        { ********** Multiply Fp27 with Zetta ***************** }
procedure _Mul_Fp27_Zetta(const Value: Fp27Int; var Result: Fp27Int);
// Value should be different than Result
begin
Result.Tower := Value.Tower;
Result.b := Value.a;
Result.c := Value.b;
_Mul_Fp9_By_V(Value.c, Result.a);
end;




        { ********** Get Square of an Fp27 ***************** }
procedure _Sqr_Fp27(const Value: Fp27Int; var Result: Fp27Int);
var tmp: array [0 .. 5] of Fp9Int;
begin
Result.Tower := Value.Tower;
  { Using Karatsuba Squering }
_Sqr_Fp9(Value.a, tmp[0]);
_Sqr_Fp9(Value.b, tmp[1]);
_Sqr_Fp9(Value.c, tmp[2]);
_Add_Fp9(Value.a, Value.b, tmp[3]);
_Add_Fp9(Value.a, Value.c, tmp[4]);
_Add_Fp9(Value.b, Value.c, Result.a);
_Sqr_Fp9(Result.a, Result.a);
_Sub_Fp9(Result.a, tmp[1], Result.a);
_Sub_Fp9(Result.a, tmp[2], Result.a);
_Mul_Fp9_By_V(Result.a, tmp[5]);
_Add_Fp9(tmp[5], tmp[0], Result.a);
_Sqr_Fp9(tmp[3], tmp[3]);
_Sub_Fp9(tmp[3], tmp[0], tmp[3]);
_Sub_Fp9(tmp[3], tmp[1], tmp[3]);
_Mul_Fp9_By_V(tmp[2], Result.b);
_Add_Fp9(Result.b, tmp[3], Result.b);
_Sqr_Fp9(tmp[4], Result.c);
_Sub_Fp9(Result.c, tmp[0], Result.c);
_Sub_Fp9(Result.c, tmp[2], Result.c);
_Add_Fp9(Result.c, tmp[1], Result.c);
end;

        { **********   Compare two Fp27 integers ************* }
function _Equals_Fp27(const Left, Right: Fp27Int): Boolean;
begin
Result := (_Equals_Fp9(Left.a, Right.a)) and (_Equals_Fp9(Left.b, Right.b))and (_Equals_Fp9(Left.c, Right.c));
end;



        { ********** Raise an Fp27 to " p " power ************ }
procedure _Pow_Fp27_P(const Value: Fp27Int; var Result: Fp27Int);
var i: integer;
begin
Result.Tower := Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.c.SetTowerParams(Value.Tower);
Result.a:=value.a.toPowerP;
_Mul_Fp_Fp9(Value.Tower^.FrobeniusP_Const[2].a,value.b.toPowerP,Result.b);
_Mul_FP9_By_V(Result.b,Result.b);
_Mul_Fp_Fp9(Value.Tower^.FrobeniusP_Const[3].a,value.c.toPowerP,Result.c);
_Mul_FP9_By_V(Result.c,Result.c);
_Mul_FP9_By_V(Result.c,Result.c);
end;

        { ********** Raise an Fp27 to "p^i " power ************ }
procedure _Pow_Fp27_P_i(const Value: Fp27Int; power: integer;var Result: Fp27Int);
var i: integer;
begin
Result:=Value;
for i := 0 to power - 1 do begin
                           Result.a:=Result.a.toPowerP;
                           _Mul_Fp_Fp9(Result.Tower^.FrobeniusP_Const[2].a,Result.b.toPowerP,Result.b);
                           _Mul_FP9_By_V(Result.b,Result.b);
                           _Mul_Fp_Fp9(Result.Tower^.FrobeniusP_Const[3].a,Result.c.toPowerP,Result.c);
                           _Mul_FP9_By_V(Result.c,Result.c);
                           _Mul_FP9_By_V(Result.c,Result.c);;
                           end;
end;

 { ********** Raise an Fp27 to " p^9 " power (Faster since a^P^9=a if a in Fp3) ************ }
procedure _Pow_Fp27_P9(const Value: Fp27Int; var Result: Fp27Int);
var i: integer;
begin
Result.Tower := Value.Tower;
Result.a.SetTowerParams(Value.Tower);
Result.b.SetTowerParams(Value.Tower);
Result.c.SetTowerParams(Value.Tower);
Result.a:=value.a;
_Mul_FP_FP9(Value.Tower.FieldParam.two_expq_1div27,Value.b ,Result.b);
_Mul_FP_FP9(Value.Tower.FieldParam.two_expq_1div27.Sqr,Value.c ,Result.c);
end;

        { ********** Get inverse of an Fp27 ***************** }
procedure _Inv_Fp27(const Value: Fp27Int; var Result: Fp27Int);
var t: array [0 .. 9] of Fp9Int;
begin
Result.Tower := Value.Tower;
_Sqr_Fp9(Value.a, t[0]);
_Sqr_Fp9(Value.b, t[1]);
_Sqr_Fp9(Value.c, t[2]);
_Mul_Fp9(Value.a, Value.b, t[3]);
_Mul_Fp9(Value.a, Value.c, t[4]);
_Mul_Fp9(Value.b, Value.c, t[5]);
_Mul_Fp9_By_V(t[5], t[7]);
_Sub_Fp9(t[0], t[7], t[7]);
_Mul_Fp9_By_V(t[2], t[8]);
_Sub_Fp9(t[8], t[3], t[8]);
_Sub_Fp9(t[1], t[4], t[9]);
_Mul_Fp9(Value.a, t[7], t[6]);
_Mul_Fp9_By_V(Value.c, t[0]);
_Mul_Fp9(t[0], t[8], t[0]);
_Add_Fp9(t[6], t[0], t[6]);
_Mul_Fp9_By_V(Value.b, t[0]);
_Mul_Fp9(t[0], t[9], t[0]);
_Add_Fp9(t[6], t[0], t[6]);
_Inv_Fp9(t[6], t[6]);
_Mul_Fp9(t[7], t[6], Result.a);
_Mul_Fp9(t[8], t[6], Result.b);
_Mul_Fp9(t[9], t[6], Result.c);
end;

        { ********** Get negative of an Fp27 ***************** }
procedure _Neg_Fp27(const Value: Fp27Int; var Result: Fp27Int);
begin
Result.Tower := Value.Tower;
_Neg_Fp9(Value.a, Result.a);
_Neg_Fp9(Value.b, Result.b);
_Neg_Fp9(Value.c, Result.c);
end;

(*        { ********** Conjugate an Fp27 ***************** }
procedure _Conjugate_Fp27(const Value: Fp27Int; var Result: Fp27Int);
begin
  Result.Tower := Value.Tower;
  _Conjugate_Fp9(Value.a, Result.a);
  _Conjugate_Fp9(Value.b, Result.b);
  _Neg_Fp9(Result.b, Result.b);
  _Conjugate_Fp9(Value.c, Result.c);

end;*)

        { ********** Raise an Fp27 to a LInt power ************ }
procedure _Pow_Fp27(const Value: Fp27Int; Exponent: LInt; var Result: Fp27Int;Poweringmode :FpPoweringMode=pmNormal);
var  i: word;
     tmp,zi: Fp27Int;
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
                                                            _Sqr_Fp27(Result, Result);
                                                            if _Is_BitSet_At(Exponent, i) then _Mul_Fp27(Result, tmp, Result);
                                                            end;
              end;

  end;
if _IsNeg(Exponent) then  _Inv_Fp27(Result,Result);
end;

        { ****** Convert an Fp27 to a Decimal String ********** }
function _Fp27_To_DecimalString(const Value: Fp27Int): string;
begin
if Value.c.IsZero then begin
    if Value.b.IsZero then begin
                           if Value.a.IsZero then
                           Result := '0'
                           else Result := Value.a.ToDecimalString;
                           end
    else begin
         if Value.b.IsOne then begin
            if Value.a.IsZero then Result := 'w'
            else Result := 'w +' + Value.a.ToDecimalString;
         end
         else  begin
              if Value.a.IsZero then Result := '(' + Value.b.ToDecimalString + ') * w'
              else Result := '(' + Value.b.ToDecimalString + ') * w +' +Value.a.ToDecimalString;
        end;
    end
  end
  else begin
    if Value.c.IsOne then begin
                          if Value.b.IsZero then begin
                          if Value.a.IsZero then Result := 'w^2'
                          else Result := 'w^2 + ' + Value.a.ToDecimalString;
                          end
      else begin
        if Value.b.IsOne then begin
                              if Value.a.IsZero then Result := 'w^2 + w'
                              else Result := 'w^2 + w +' + Value.a.ToDecimalString;
                              end
        else begin
             if Value.a.IsZero then Result := 'w^2 +(' + Value.b.ToDecimalString + ') * w'
             else  Result := 'w^2 + (' + Value.b.ToDecimalString + ') * w +' +Value.a.ToDecimalString;
             end;
        end
    end
    else begin
      if Value.b.IsZero then begin
                             if Value.a.IsZero then Result := '(' + Value.c.ToDecimalString + ') * w^2'
                             else Result := '(' + Value.c.ToDecimalString + ') * w^2 + ' +Value.a.ToDecimalString;
                             end
      else begin
           if Value.b.IsOne then begin
           if Value.a.IsZero then Result := '(' + Value.c.ToDecimalString + ') * w^2 +' + 'w'
           else Result := '(' + Value.c.ToDecimalString + ') * w^2 +' + ' w +' +Value.a.ToDecimalString;
        end
        else begin
             if Value.a.IsZero then Result := '(' + Value.c.ToDecimalString + ') * w^2 +' + '(' +Value.b.ToDecimalString + ') * w'
             else Result := '(' + Value.c.ToDecimalString + ') * w^2 +' + '(' +Value.b.ToDecimalString + ') * w +' + Value.a.ToDecimalString;
             end;
          end
      end;
end;
end;

        { ****** Convert an Fp27 to a Hexadecimal String ******* }
function _Fp27_To_HexString(const Value: Fp27Int): string;
begin
  if Value.c.IsZero then begin
                         if Value.b.IsZero then begin
                                                if Value.a.IsZero then Result := '0'
                                                else Result := Value.a.ToHexString;
                                                end
                         else begin
                              if Value.b.IsOne then begin
                                                    if Value.a.IsZero then Result := 'w'
                                                    else Result := 'w +' + Value.a.ToHexString;
                                                    end
                              else begin
                                   if Value.a.IsZero then Result := '(' + Value.b.ToHexString + ') * w'
                                   else Result := '(' + Value.b.ToHexString + ') * w +' + Value.a.ToHexString;
                                   end;
                              end
                         end
  else  begin
        if Value.c.IsOne then begin
                              if Value.b.IsZero then begin
                                                     if Value.a.IsZero then Result := 'w^2'
                                                     else Result := 'w^2 + ' + Value.a.ToHexString;
                                                     end
                              else begin
                                   if Value.b.IsOne then begin
                                                         if Value.a.IsZero then Result := 'w^2 + w'
                                                         else Result := 'w^2 + w +' + Value.a.ToHexString;
                                                         end
                                   else begin
                                        if Value.a.IsZero then Result := 'w^2 +(' + Value.b.ToHexString + ') * w'
                                        else Result := 'w^2 + (' + Value.b.ToHexString + ') * w +' +Value.a.ToHexString;
                                        end;
                                   end
                              end
        else begin
             if Value.b.IsZero then begin
                                    if Value.a.IsZero then Result := '(' + Value.c.ToHexString + ') * w^2'
                                    else Result := '(' + Value.c.ToHexString + ') * w^2 + ' +Value.a.ToHexString;
                                    end
             else begin
                  if Value.b.IsOne then begin
                                        if Value.a.IsZero then Result := '(' + Value.c.ToHexString + ') * w^2 +' + 'w'
                                        else Result := '(' + Value.c.ToHexString + ') * w^2 +' + ' w +' +
                                        Value.a.ToHexString;
                                        end
                  else begin
                       if Value.a.IsZero then Result := '(' + Value.c.ToHexString + ') * w^2 +' + '(' +Value.b.ToHexString + ') * w'
                       else Result := '(' + Value.c.ToHexString + ') * w^2 +' + '(' +Value.b.ToHexString + ') * w +' + Value.a.ToHexString;
                       end;
                  end
             end;
        end;
end;

{ ****** Convert String to an Fp27 (Decimal/Hex)******* }
procedure _Fp27_From_Strings(const valueA, valueB, valueC, valueD, valueE,
  ValueF, ValueG, ValueH, ValueI: String;
  var Result: Fp27Int);
begin
  Result.a.SetFromStrings(valueA,ValueB,ValueC);
  Result.b.SetFromStrings(valueD,ValueE,ValueF);
  Result.c.SetFromStrings(ValueG,ValueH,ValueI);
end;

{ ******************************************************************************* }
/// Definitions of an Fp27 integer operators and functions
{ ******************************************************************************* }

class operator Fp27Int.Add(const Left, Right: Fp27Int): Fp27Int;
begin
  _Add_Fp27(Left, Right, Result);
end;

class operator Fp27Int.Subtract(const Left, Right: Fp27Int): Fp27Int;
begin
  _Sub_Fp27(Left, Right, Result);
end;

class operator Fp27Int.Multiply(const Left, Right: Fp27Int): Fp27Int;
begin
  _Mul_Fp27(Left, Right, Result);
end;

class operator Fp27Int.Multiply(const Left: Fp9Int; Right: Fp27Int): Fp27Int;
begin
  _Mul_Fp9_Fp27(Left, Right, Result);
end;

function Fp27Int.MultiplyByZetta: Fp27Int;
begin
  _Mul_Fp27_Zetta(Self, Result);
end;

function Fp27Int.Sqr: Fp27Int;
begin
  _Sqr_Fp27(Self, Result);
end;

class operator Fp27Int.Equal(const Left, Right: Fp27Int): Boolean;
begin
  Result := _Equals_Fp27(Left, Right);
end;

function Fp27Int.Inverse: Fp27Int;
begin
  _Inv_Fp27(Self, Result);
end;

function Fp27Int.IsOne: Boolean;
begin
  Result := (a.IsOne) and (b.IsZero) and (c.IsZero);
end;

function Fp27Int.IsZero: Boolean;
begin
  Result := (a.IsZero) and (b.IsZero) and (c.IsZero);
end;

class operator Fp27Int.Negative(const Value: Fp27Int): Fp27Int;
begin
  _Neg_Fp27(Value, Result);
end;

class operator Fp27Int.NotEqual(const Left, Right: Fp27Int): Boolean;
begin
  Result := not _Equals_Fp27(Left, Right)
end;

function Fp27Int.Pow(Exponent: LInt;Poweringmode:FpPoweringMode=pmNormal): Fp27Int;
begin
  _Pow_Fp27(Self, Exponent, Result,Poweringmode);
end;

procedure Fp27Int.SetTowerParams(TParam: PtrTowerParams18);
begin
  Tower := TParam;
  a.SetTowerParams(TParam);
  b.SetTowerParams(TParam);
  c.SetTowerParams(TParam);
end;

procedure Fp27Int.SetToZero;
begin
  a.SetToOne;
  b.SetToZero;
  c.SetToZero;
end;

procedure Fp27Int.SetFromStrings(a, b, c, d, e, f, g, h, i, j, k, l: String);
begin
  _Fp27_From_Strings(a, b, c, d, e, f, g, h, i, Self);
end;

procedure Fp27Int.SetToOne;
begin
  a.SetToOne;
  b.SetToZero;
  c.SetToZero;
end;

procedure Fp27Int.SetToRandom;
begin
  a.SetToRandom;
  b.SetToRandom;
  c.SetToRandom;
end;

function Fp27Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.a.BitLength mod 8) * 8;
Setlength(Result,size*27);
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

Move(c.a.a.Data.i8[0],Result[18*size],size);
Move(c.a.b.Data.i8[0],Result[19*size],size);
Move(c.a.c.Data.i8[0],Result[20*size],size);
Move(c.b.a.Data.i8[0],Result[21*size],size);
Move(c.b.b.Data.i8[0],Result[22*size],size);
Move(c.b.c.Data.i8[0],Result[23*size],size);
Move(c.c.a.Data.i8[0],Result[24*size],size);
Move(c.c.b.Data.i8[0],Result[25*size],size);
Move(c.c.c.Data.i8[0],Result[26*size],size);
end;

function Fp27Int.ToDecimalString: string;
begin
  Result := _Fp27_To_DecimalString(Self)
end;

function Fp27Int.ToHexString: string;
begin
  Result := _Fp27_To_HexString(Self)
end;

end.
