unit Fp3Arithmetic;

interface
uses System.SysUtils,VCL.dialogs, VLargeIntegers,LargeIntegers;

type

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp3. Fp3 is the cubic extention of Fp with respect to the
   irreducible polynômial u^3+Beta (u^3=Beta). Elements are on the form a+u*b+u^2*c (a,b,c from Fp).

********************************************************************************************}


             { Definition of the Fp3 field Extension of the field FP modulo the irriductibe polynomial X^2+Beta=0 }

  Fp3Int=record
  Field:PtrFieldParams;
  a,b,c:LInt;
    public
    function toHexString: string;
    function toDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse:Fp3Int;
    function Sqr:Fp3Int;
    function Sqrt:Fp3Int;
    function toPowerP:Fp3Int;
    function Pow(Exponent: LInt): Fp3Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function IsASquare:boolean;
    procedure SetToRandom;
    procedure SetToZero;
    procedure SetToOne;
    procedure SetFieldParams(FieldParams:PtrFieldParams);
    procedure SetFormString(const Value: string);
    class operator Add(const Left, Right: Fp3Int): Fp3Int;
    class operator Add(const Left:LInt; Right: Fp3Int): Fp3Int;
    class operator Add(const Left: Fp3Int;Right:LInt): Fp3Int;
    class operator Subtract(const Left, Right: Fp3Int): Fp3Int;
    class operator Subtract(const Left:LInt; Right: Fp3Int): Fp3Int;
    class operator Subtract(const Left:Fp3Int; Right: LInt): Fp3Int;
    class operator Multiply(const Left, Right: Fp3Int): Fp3Int;
    class operator Multiply(const Left: Fp3Int; Right: Word): Fp3Int;
    class operator Multiply(Left: Word; const Right: Fp3Int): Fp3Int;
    class operator Multiply(Left: LInt; const Right: Fp3Int): Fp3Int;
    class operator Modulus(const Left: Fp3Int; Right: LInt): Fp3Int;
    class operator Negative(const Value: Fp3Int): Fp3Int;
    class operator Equal(const Left, Right: Fp3Int): Boolean;
    class operator Equal(const Left:Fp3Int;Right:LInt): Boolean;
    class operator NotEqual(const Left, Right: Fp3Int): Boolean;
    class operator GreaterThan(const Left, Right: Fp3Int): Boolean;
    class operator GreaterThanOrEqual(const Left, Right: Fp3Int): Boolean;
    class operator LessThan(const Left, Right: Fp3Int): Boolean;
    class operator LessThanOrEqual(const Left, Right: Fp3Int): Boolean;
  end;

  procedure _Add_Fp3(Left,Right:Fp3Int; var Result:Fp3Int);
  procedure _Add_FP_Fp3(Left:LInt; Right: Fp3Int; var Result:Fp3Int);
  procedure _Add_Fp3_FP(Left:Fp3Int; Right: LInt; var Result:Fp3Int);
  procedure _Sub_Fp3(Left,Right:Fp3Int;var Result: Fp3Int);
  procedure _Sub_Fp3_FP(Left:Fp3Int; Right: LInt; var Result:Fp3Int);
  function  _Compare_Fp3(const Left,Right:Fp3Int):shortint;
  function  _Compare_Fp3_FP(const Left:Fp3Int;Right:LInt):boolean;
  procedure _Mul_Fp3(const Left,Right:Fp3Int;var Result:Fp3Int);
  procedure _Mul_Fp3_FP(const Left: Fp3Int; Right: LInt; var Result: Fp3Int);
  procedure _Mul_FP_Fp3(const Left: LInt; Right: Fp3Int; var Result: Fp3Int);
  procedure _Mul_Fp3_Word(const Left: Fp3Int; Right: Word; var Result: Fp3Int);
  procedure _Mul_Word_Fp3(const Left: Word; Right: Fp3Int; var Result: Fp3Int);
  procedure _Neg_Fp3(const Value: Fp3Int; var Result: Fp3Int);
  procedure _Mod_Fp3_FP(const Left: Fp3Int; Right: LInt; var Result: Fp3Int);
  procedure _Sqr_Fp3(const Value:Fp3Int;var Result:Fp3Int);
  procedure _Inv_Fp3(const value:Fp3Int;var Result:Fp3Int);
  procedure _Sqrt_Fp3(const Value:Fp3Int;var Result:Fp3Int);
  procedure _Pow_Fp3(const value:Fp3Int;Exponent:LInt;var Result:Fp3Int);
  function  _Fp3_To_HexString(const value:Fp3Int):string;
  procedure _Fp3_From_String(const value:String;var Result:Fp3Int);
  function  _Is_Fp3_ASquare(const value:Fp3Int):boolean;
  function  _Is_Fp3_Null(const value:Fp3Int):boolean;
  function  InitFieldParamsFp3(Beta:integer;P:Lint):PtrFieldParams;
  procedure _NoMod_Add_Fp3(Left,Right:Fp3Int; var Result:Fp3Int);
  procedure _NoMod_Sub_Fp3(Left,Right:Fp3Int;var Result: Fp3Int);
  procedure _PowP_Fp3(Value: Fp3Int; var Result:Fp3Int);
  procedure _Mul_Fp3_u(const value: Fp3Int; var Result: Fp3Int);
  procedure _Div_Fp3_u(const value: Fp3Int; var Result: Fp3Int);

implementation

function InitFieldParamsFp3(Beta:integer;P:Lint):PtrFieldParams;
var tmpsru,tmpsru2:HLInt;
begin
New(Result);
Result.Beta:=Beta;
Result^.p:=p;
Result^.p_1_div_2:=(p-1) shr 1;
Result^.inv_2:=LInt(2).InversModulo(p);
New(Result.MontgomeryData);
InitMontgomeryStruct(P,Result.MontgomeryData);
//Result^.sru:=Lint(Beta).PowMod((p-1)/6,p,Result.MontgomeryData); // For Square Root Computation
Result.sru:=Lint(Beta).PowMod((p-1)/6,p,Result.MontgomeryData); // For Square Root Computation
_VHCopy2_LInt(Result.sru,tmpsru);
VLargeIntegers._Sqr_LInt(tmpsru,tmpsru2);
VLargeIntegers._Mod_LInt(tmpsru2,p,tmpsru2);
_VHCopy_LInt(tmpsru2,Result.sru2);
_VHCopy_LInt(tmpsru,Result.sru);
Result.pmod8:=p mod 8;
Result.p_5_div_8:=(p-5)/8;
Result.p_3_div_4:=(p-3)/4;
_VHCopy2_LInt(Result.p,tmpsru);
tmpsru:=tmpsru.Pow(9);
tmpsru2:=(tmpsru-1)/27;
tmpsru:=(tmpsru-5)/8;
tmpsru:=HLint(2).PowMod(tmpsru,Result.p);
tmpsru2:=HLint(2).PowMod(tmpsru2,Result.p);
_VHCopy_LInt(tmpsru,result.two_expq_div8);
_VHCopy_LInt(tmpsru2,result.two_expq_1div27);
end;


{*******************************************************************************}
            ///      Procedures for Fp3 Arithmetic
{*******************************************************************************}

        {**********   Add two Fp3 integers *****************}
        // Without Modular reduction (for optimization in some cases)
procedure _NoMod_Add_Fp3(Left,Right:Fp3Int; var Result:Fp3Int);
begin
Result.Field:=Left.Field;
_Add_LInt(Left.a,Right.a,Result.a);
_Add_LInt(Left.b,Right.b,Result.b);
_Add_LInt(Left.c,Right.c,Result.c);
end;

        {**********   Add two Fp3 integers *****************}
procedure _Add_Fp3(Left,Right:Fp3Int; var Result:Fp3Int);
begin
Result.Field:=Left.Field;
_Add_LInt(Left.a,Right.a,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Add_LInt(Left.b,Right.b,Result.b);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
_Add_LInt(Left.c,Right.c,Result.c);
_Mod_LInt(Result.c,Left.Field^.p,Result.c);
end;
        {**********   Add LInt with Fp3 integer***********}
procedure _Add_FP_Fp3(Left:LInt; Right: Fp3Int; var Result:Fp3Int);
begin
Result.Field:=Right.Field;
_Add_LInt(Right.a,Left,Result.a);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
Result.b:=Right.b;
Result.c:=Right.c;
end;
        {**********   Add Fp3 integer with LInt **********}
procedure _Add_Fp3_FP(Left:Fp3Int; Right: LInt; var Result:Fp3Int);
begin
Result.Field:=Left.Field;
_Add_LInt(Left.a,Right,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
Result.b:=Left.b;
Result.c:=Left.c;
end;
        {**********   Substract two Fp3 integers ***********}
        // Without Modular reduction (for optimization in some cases)
procedure _NoMod_Sub_Fp3(Left,Right:Fp3Int;var Result: Fp3Int);
begin
Result.Field:=Left.Field;
_Sub_LInt(Left.a,Right.a,Result.a);
_Sub_LInt(Left.b,Right.b,Result.b);
_Sub_LInt(Left.c,Right.c,Result.c);
end;

        {**********   Substract two Fp3 integers ***********}
procedure _Sub_Fp3(Left,Right:Fp3Int;var Result: Fp3Int);
begin
Result.Field:=Left.Field;
_Sub_LInt(Left.a,Right.a,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Sub_LInt(Left.b,Right.b,Result.b);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
_Sub_LInt(Left.c,Right.c,Result.c);
_Mod_LInt(Result.c,Left.Field^.p,Result.c);
end;
        {*********   Substract LInt with Fp3 integer *****}
procedure _Sub_FP_Fp3(const Left:LInt; Right: Fp3Int; var Result:Fp3Int);
begin
Result.Field:=Right.Field;
Result.b:=Right.b;
Result.c:=Right.c;
_Sub_LInt(Left,Right.a,Result.a);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
end;
        {******** Substract Fp3 integer with LInt ********}
procedure _Sub_Fp3_FP(Left:Fp3Int; Right: LInt; var Result:Fp3Int);
begin
Result.Field:=Left.Field;
Result.b:=Left.b;
Result.c:=Left.c;
_Sub_LInt(Left.a,Right,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
end;

        {**********   Compare two Fp3 integers *************}
function _Compare_Fp3(const Left,Right:Fp3Int):shortint;
begin
if (_Compare_LInt(Left.a,right.a)=0) and(_Compare_LInt(left.b,right.b)=0) and(_Compare_LInt(left.c,right.c)=0) then result:=0
else result:=-1;
end;
        {******** Compare Fp3 integer with LInt **********}
function _Compare_Fp3_FP(const Left:Fp3Int;Right:LInt):boolean;
begin
Result:=(Left.b.Data.i16[-2]=0)and(Left.c.Data.i16[-2]=0)and(_Compare_LInt(Left.a,Right)=0);
end;

        {********** Multiply two Fp3 integers *****************}
procedure _Mul_Fp3(const Left,Right:Fp3Int;var Result:Fp3Int);
var tmp:array[0..10] of LInt;
begin
Result.Field:=Left.Field;
        { Multiplication Using Karatsuba Multiplication  }
Result.Field:=Left.Field;
_Mul_LInt(Left.a,Right.a,tmp[0]);
_Mul_LInt(Left.b,Right.b,tmp[1]);
_Mul_LInt(Left.c,Right.c,tmp[2]);
_Mod_LInt(tmp[0],Left.Field.p,tmp[0]);
_Mod_LInt(tmp[1],Left.Field.p,tmp[1]);
_Mod_LInt(tmp[2],Left.Field.p,tmp[2]);
_Add_LInt(Left.b,Left.c,tmp[3]);
_Add_LInt(Right.b,Right.c,tmp[4]);
_Add_LInt(Left.a,Left.b,tmp[5]);
_Add_LInt(Right.a,Right.b,tmp[6]);
_Add_LInt(Left.a,Left.c,tmp[7]);
_Add_LInt(Right.a,Right.c,tmp[8]);
_Mul_LInt(tmp[3],tmp[4],Result.a);
_Sub_LInt(Result.a,tmp[1],Result.a);
_Sub_LInt(Result.a,tmp[2],Result.a);
_Mul_LInt(Result.a,Lint(Left.Field.Beta),tmp[9]);
_Add_LInt(tmp[0],tmp[9],Result.a);
_Mod_LInt(Result.a,Left.Field.p,Result.a);

_Mul_LInt(tmp[5],tmp[6],Result.b);
_Sub_LInt(Result.b,tmp[0],Result.b);
_Sub_LInt(Result.b,tmp[1],Result.b);
_Mul_LInt(tmp[2],Lint(Left.Field.Beta),tmp[9]);
_Add_LInt(Result.b,tmp[9],Result.b);
_Mod_LInt(Result.b,Left.Field.p,Result.b);

_Mul_LInt(tmp[7],tmp[8],Result.c);
_Sub_LInt(Result.c,tmp[0],Result.c);
_Sub_LInt(Result.c,tmp[2],Result.c);
_Add_LInt(Result.c,tmp[1],Result.c);
_Mod_LInt(Result.c,Left.Field.p,Result.c);
end;


        {********** Multiply Fp3 with LInt *****************}
procedure _Mul_Fp3_FP(const Left: Fp3Int; Right: LInt; var Result: Fp3Int);
var x,y,z:Lint;
begin
Result.Field:=Left.Field;
_Mul_LInt(Left.a,Right,x);
_Mul_LInt(Left.b,Right,y);
_Mul_LInt(Left.c,Right,z);
_Mod_LInt(x,Left.Field^.p,Result.a);
_Mod_LInt(y,Left.Field^.p,Result.b);
_Mod_LInt(z,Left.Field^.p,Result.c);
end;
        {********** Multiply LInt with Fp3 *****************}
procedure _Mul_FP_Fp3(const Left: LInt; Right: Fp3Int; var Result: Fp3Int);
var x,y,z:lint;
begin
Result.Field:=Right.Field;
_Mul_LInt(Left,Right.a,x);
_Mul_LInt(Left,Right.b,y);
_Mul_LInt(Left,Right.c,z);
_Mod_LInt(x,Right.Field^.p,Result.a);
_Mod_LInt(y,Right.Field^.p,Result.b);
_Mod_LInt(z,Right.Field^.p,Result.c);
end;
        {********** Multiply Fp3 with Word *****************}
procedure _Mul_Fp3_Word(const Left: Fp3Int; Right: Word; var Result: Fp3Int);
begin
Result.Field:=Left.Field;
_Mul_LInt(Left.a,LInt(Right),Result.a);
_Mul_LInt(Left.b,LInt(Right),Result.b);
_Mul_LInt(Left.c,LInt(Right),Result.c);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
_Mod_LInt(Result.c,Left.Field^.p,Result.c);
end;
        {********** Multiply Word with Fp3 *****************}
procedure _Mul_Word_Fp3(const Left: Word; Right: Fp3Int; var Result: Fp3Int);
begin
Result.Field:=Right.Field;
_Mul_LInt(LInt(Left),Right.a,Result.a);
_Mul_LInt(LInt(Left),Right.b,Result.b);
_Mul_LInt(LInt(Left),Right.c,Result.c);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
_Mod_LInt(Result.b,Right.Field^.p,Result.b);
_Mod_LInt(Result.c,Right.Field^.p,Result.c);
end;

        {********** Multiply Fp3 by u *****************} //used for FP6 Sqrt
procedure _Mul_Fp3_u(const value: Fp3Int; var Result: Fp3Int);
var x:Lint;
begin
Result.Field:=Value.Field;
_HCopy_LInt(value.b,x);
_HCopy_LInt(value.a,result.b);
_Mul_LInt(value.c,Lint(value.Field.Beta),Result.a);
_Mod_LInt(Result.a,value.Field.p, Result.a);
_HCopy_LInt(x,result.c);
end;

        {********** Divide Fp3 by u *****************}  //used for FP6 Sqrt
procedure _Div_Fp3_u(const value: Fp3Int; var Result: Fp3Int);
var a,x:lint;
begin
Result.Field:=Value.Field;
_HCopy_LInt(value.a,x);
_HCopy_LInt(value.b,result.a);
_HCopy_LInt(value.c,result.b);
a.Data.i16[-2]:=1;
if value.Field.Beta<0 then a.Data.i16[-1]:=$ffff
else a.Data.i16[-1]:=0;
a.Data.i32[0]:=abs(value.Field.Beta);
_Mul_LInt(x,a.InversModulo(value.Field.p),Result.c);
end;

        {********** Get negative of an Fp3 *****************}
procedure _Neg_Fp3(const Value: Fp3Int; var Result: Fp3Int);
begin
result.Field:=Value.Field;
if Value.a.Data.i16[-2]<>0 then  _Sub_LInt(Value.Field^.p,Value.a,Result.a) else Result.a.Data.i16[-2]:=0;
if Value.b.Data.i16[-2]<>0 then  _Sub_LInt(Value.Field^.p,Value.b,Result.b)else Result.b.Data.i16[-2]:=0;
if Value.c.Data.i16[-2]<>0 then  _Sub_LInt(Value.Field^.p,Value.c,Result.c)else Result.c.Data.i16[-2]:=0;
end;
        {********** Get modulo of an Fp3 *****************}
procedure _Mod_Fp3_FP(const Left: Fp3Int; Right: LInt; var Result: Fp3Int);
begin
Result.Field:=Left.Field;
_Mod_LInt(Left.a,Right,Result.a);
_Mod_LInt(Left.b,Right,Result.b);
_Mod_LInt(Left.c,Right,Result.c);
end;

        {********** Get Square of an Fp3 *****************}
procedure _Sqr_Fp3(const Value:Fp3Int;var Result:Fp3Int);
var tmp:array[0..5]of LInt;
begin
Result.Field:=Value.Field;
      {       Using Karatsuba Squering    }
_Sqr_LInt(Value.a,tmp[0]);
_Sqr_LInt(Value.b,tmp[1]);
_Sqr_LInt(Value.c,tmp[2]);
_Add_LInt(Value.a,Value.b,tmp[3]);
_Add_LInt(Value.a,Value.c,tmp[4]);
_Add_LInt(Value.b,Value.c,Result.a);
_Sqr_LInt(Result.a,tmp[5]);
_Sub_LInt(tmp[5],tmp[1],Result.a);
_Sub_LInt(Result.a,tmp[2],Result.a);
_Mul_LInt(Lint(value.Field.Beta),Result.a,tmp[5]);
_Add_LInt(tmp[5],tmp[0],Result.a);
_Mod_LInt(Result.a,Value.Field.p,Result.a);
_Sqr_LInt(tmp[3],tmp[5]);
_Sub_LInt(tmp[5],tmp[0],tmp[3]);
_Sub_LInt(tmp[3],tmp[1],tmp[3]);
_Mul_LInt(Lint(value.Field.Beta),tmp[2],Result.b);
_Add_LInt(Result.b,tmp[3],Result.b);
_Mod_LInt(Result.b,Value.Field.p,Result.b);
_Sqr_LInt(tmp[4],Result.c);
_Sub_LInt(Result.c,tmp[0],Result.c);
_Sub_LInt(Result.c,tmp[2],Result.c);
_Add_LInt(Result.c,tmp[1],Result.c);
_Mod_LInt(Result.c,Value.Field.p,Result.c);
end;

        {********** Get inverse of an Fp3 *****************}
procedure _Inv_Fp3(const value:Fp3Int;var Result:Fp3Int);
var t:array[0..10] of LInt;
begin
Result.Field:=Value.Field;
_Sqr_LInt(value.a,t[0]);
_Sqr_LInt(value.b,t[1]);
_Sqr_LInt(value.c,t[2]);
_Mul_LInt(value.a,value.b,t[3]);
_Mul_LInt(value.a,value.c,t[4]);
_Mul_LInt(value.b,value.c,t[5]);
_Mul_LInt(Lint(value.Field.Beta),t[5],t[7]);
_Sub_LInt(t[0],t[7],t[7]);
_Mul_LInt(Lint(value.Field.Beta),t[2],t[8]);
_Sub_LInt(t[8],t[3],t[8]);
_Sub_LInt(t[1],t[4],t[9]);
_Mul_LInt(value.a,t[7],t[6]);
_Mul_LInt(Lint(value.Field.Beta),value.c,t[0]);
_Mul_LInt(t[0],t[8],t[10]);
_Add_LInt(t[6],t[10],t[6]);
_Mul_LInt(Lint(value.Field.Beta),value.b,t[0]);
_Mul_LInt(t[0],t[9],t[10]);
_Add_LInt(t[6],t[10],t[6]);
_Inv_Mod_LInt(t[6],Value.Field.p,t[6]);
_Mul_LInt(t[7],t[6],Result.a);
_Mul_LInt(t[8],t[6],Result.b);
_Mul_LInt(t[9],t[6],Result.c);
_Mod_LInt(Result.a,Value.Field.p,Result.a);
_Mod_LInt(Result.b,Value.Field.p,Result.b);
_Mod_LInt(Result.c,Value.Field.p,Result.c);
end;

    {********** Compute power to P (the modulo) of an Fp3 ************}
    // it was simply the conjugate for Fp2 !!!!!
procedure _PowP_Fp3(Value: Fp3Int; var Result:Fp3Int);
var tmp:Lint;
begin
Result.Field:=Value.Field;
if @Value<>@Result then _HCopy_LInt(Value.a,Result.a);
_Mul_LInt(Value.b,Value.Field.sru2,Result.b);
_Mod_LInt(Result.b,Value.Field.p,Result.b);
_Mul_LInt(Value.c,Value.Field.sru2,Result.c);
_Mod_LInt(Result.c,Value.Field.p,Result.c);
_Mul_LInt(Result.c,Value.Field.sru2,tmp);
_Mod_LInt(tmp,Value.Field.p,Result.c);
end;

        {********** Get Square root of an Fp3 ************}
procedure _Sqrt_Fp3(const Value:Fp3Int;var Result:Fp3Int);
var t:array[0..3] of Fp3Int;
begin
Result.Field:=Value.Field;
if not _Is_Fp3_ASquare(Value) then  raise ERangeError.Create( 'This FP3 element is not a Square .....')
else begin
     if Value.Field.pmod8=5 then begin
                                 _Add_Fp3(Value,Value,t[0]);
                                 _PowP_Fp3(t[0],t[0]);
                                 _Sqr_Fp3(t[0],t[1]);
                                 _Mul_Fp3(t[1],t[0],t[2]);
                                 _Mul_Fp3(t[1],t[2],t[1]);
                                 _PowP_Fp3(t[0],t[0]);
                                 _Add_Fp3(Value,Value,t[3]);
                                 _Mul_Fp3(t[3],t[1],t[3]);
                                 _Mul_Fp3(t[3],t[0],t[0]);
                                 _Pow_Fp3(t[0],Value.Field.p_5_div_8,t[3]);
                                 _Mul_Fp3(t[2],t[3],t[0]);
                                 _Mul_Fp3(t[0],t[0],t[1]);
                                 _Mul_Fp3(Value,t[1],t[1]);
                                 _Add_Fp3(t[1],t[1],t[1]);
                                 _Mul_Fp3(t[0],Value,t[0]);
                                 _Dec_LInt(t[1].a,1);
                                 if _IsNeg(t[1].a) then _Add_Lint(t[1].a,Value.Field.p,t[1].a);
                                 _Mul_Fp3(t[0],t[1],Result);
                                 end
     else if Value.Field.pmod8=7 then begin   {p mod4 =3}
                                      _PowP_Fp3(Value,t[0]);
                                      _Sqr_Fp3(t[0],t[1]);
                                      _Mul_Fp3(t[0],t[1],t[2]);
                                      _PowP_Fp3(t[0],t[0]);
                                      _Mul_Fp3(t[2],Value,t[3]);
                                      _Mul_Fp3(t[3],t[0],t[0]);
                                      _Pow_Fp3(t[0],Value.Field.p_3_div_4,t[2]);
                                      _Mul_Fp3(Value,t[1],t[0]);
                                      _Mul_Fp3(t[2],t[0],Result);
                                      end
     else raise ERangeError.Create( 'Can''t Compte  Square root for this element, p mod 8 is neither 5 or 7 .....')
     end;
end;
        {********** Raise an Fp3 to a LInt power ************}
procedure _Pow_Fp3(const value:Fp3Int;Exponent:LInt;var Result:Fp3Int);
var i:word;
begin
Result.Field:= Value.Field;
if (Value.b.Data.i16[-2]=0)and(Value.c.Data.i16[-2]=0) then  begin
    _Mg_Mod_Pow_LInt(Value.a,Exponent,Value.Field^.p,Value.Field.MontgomeryData,Result.a);
    Result.b.Data.i16[-2]:=0;
    Result.c.Data.i16[-2]:=0;
    end
else begin
     Result.a:=Value.a;
     Result.b:=Value.b;
     Result.c:=Value.c;
     if (Exponent>1) then
     for i:=Exponent.BitLength-2 downto 0 do begin
                                              _Sqr_Fp3(Result,Result);
                                              if _Is_BitSet_At(Exponent,i) then
                                              _Mul_Fp3(Result,Value,Result);
                                              end;
     end;
end;

{****** Convert an Fp3 to a Decimal String **********}
function _Fp3_To_DecimalString(const value:Fp3Int):string;
begin
if Value.c=0 then begin
                 if Value.b.IsZero then begin
                                  if Value.a.IsZero then Result:='0'
                                  else Result:=Value.a.ToDecimalString;
                                  end
                 else begin
                      if Value.b=1 then begin
                                      if Value.a.IsZero then Result:='u'
                                      else Result:='u +'+Value.a.ToDecimalString;
                                      end
                       else begin
                            if Value.a.IsZero then Result:='('+Value.b.ToDecimalString+') * u'
                            else Result:='('+Value.b.ToDecimalString+') * u +'+Value.a.ToDecimalString;
                            end;
                      end
                 end
else begin
     if Value.c=0 then begin
                     if Value.b.IsZero then begin
                                      if Value.a.IsZero then Result:='u^2'
                                      else Result:='u^2 + '+Value.a.ToDecimalString;
                                      end
                     else begin
                          if Value.b=1 then begin
                                          if Value.a.IsZero then Result:='u^2 + u'
                                          else Result:='u^2 + u +'+Value.a.ToDecimalString;
                                          end
                          else begin
                               if Value.a.IsZero then Result:='u^2 +('+Value.b.ToDecimalString+') * u'
                               else Result:='u^2 + ('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                               end;
                         end
                    end
    else begin
         if Value.b=0 then begin
                          if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * u^2'
                                      else Result:='('+Value.c.ToDecimalString+') * u^2 + '+Value.a.ToDecimalString;
                                      end
                          else begin
                               if Value.b=1 then begin
                                               if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * u^2 +'+ 'u'
                                               else Result:='('+Value.c.ToDecimalString+') * u^2 +'+' u +'+Value.a.ToDecimalString;
                                               end
                               else begin
                                    if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * u^2 +'+'('+Value.b.ToDecimalString+') * u'
                                    else Result:='('+Value.c.ToDecimalString+') * u^2 +'+'('+Value.b.ToDecimalString+') * u +'+Value.a.ToDecimalString;
                                    end;
                              end
                          end;
        end;
end;

        {****** Convert an Fp3 to a Hexadecimal String *******}
function _Fp3_To_HexString(const value:Fp3Int):string;
begin
if Value.c=0 then begin
                 if Value.b.IsZero then begin
                                  if Value.a.IsZero then Result:='0'
                                  else Result:=Value.a.ToHexString;
                                  end
                 else begin
                      if Value.b=1 then begin
                                      if Value.a.IsZero then Result:='u'
                                      else Result:='u +'+Value.a.ToHexString;
                                      end
                       else begin
                            if Value.a.IsZero then Result:='('+Value.b.ToHexString+') * u'
                            else Result:='('+Value.b.ToHexString+') * u +'+Value.a.ToHexString;
                            end;
                      end
                 end
else begin
     if Value.c=0 then begin
                     if Value.b.IsZero then begin
                                      if Value.a.IsZero then Result:='u^2'
                                      else Result:='u^2 + '+Value.a.ToHexString;
                                      end
                     else begin
                          if Value.b=1 then begin
                                          if Value.a.IsZero then Result:='u^2 + u'
                                          else Result:='u^2 + u +'+Value.a.ToHexString;
                                          end
                          else begin
                               if Value.a.IsZero then Result:='u^2 +('+Value.b.ToHexString+') * u'
                               else Result:='u^2 + ('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                               end;
                         end
                    end
    else begin
         if Value.b=0 then begin
                          if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * u^2'
                                      else Result:='('+Value.c.ToHexString+') * u^2 + '+Value.a.ToHexString;
                                      end
                          else begin
                               if Value.b=1 then begin
                                               if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * u^2 +'+ 'u'
                                               else Result:='('+Value.c.ToHexString+') * u^2 +'+' u +'+Value.a.ToHexString;
                                               end
                               else begin
                                    if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * u^2 +'+'('+Value.b.ToHexString+') * u'
                                    else Result:='('+Value.c.ToHexString+') * u^2 +'+'('+Value.b.ToHexString+') * u +'+Value.a.ToHexString;
                                    end;
                              end
                          end;
        end;
end;
        {****** Convert String to an Fp3 (Decimal/Hex)*******}
procedure _Fp3_From_String(const value:String;var Result:Fp3Int);
var i:integer;
    s1,s2,s3,Val:string;
    valid:boolean;
begin
Val:=value;
if pos('+u^2*',Val)=0 then s3:=''
else begin
     Valid:=true;
     for i:=pos('+u^2*',Val)+5 to length(Val) do if not(Val[i]in ['-','x','$','0'..'9','a'..'f','A'..'F']) then valid:=false;
     if not valid then raise Exception.Create('Valeur Invalide pout un type Fp3Int')
     else begin
          S3:=Copy(Val,pos('+u^2*',Val)+5,length(Val));
          Val:=Copy(Val,1,pos('+u^2*',Val)-1);
          end;
     end;
if pos('+u*',Val)=0 then s2:=''
else begin
     Valid:=true;
     for i:=pos('+u*',Val)+3 to length(Val) do if not(Val[i]in ['-','x','$','0'..'9','a'..'f','A'..'F']) then valid:=false;
     if not valid then raise Exception.Create('Valeur Invalide pout un type Fp3Int')
     else begin
          S2:=Copy(Val,pos('+u*',Val)+3,length(Val));
          Val:=Copy(Val,1,pos('+u*',Val)-1);
          end;
     end;
s1:=Val;
if S3='' then Result.c:=0 Else Result.c:=s3;
if S2='' then Result.b:=0 Else Result.b:=s2;
if Val='' then Result.a:=0 Else Result.a:=Val;
end;

       {****** Test if an Fp3 is a Sequare *******}
function _Is_Fp3_ASquare(const value:Fp3Int):boolean;
var t:array[0..1] of Fp3Int;
begin
_PowP_Fp3(Value,t[0]);
_Mul_Fp3(Value,t[0],t[1]);
_PowP_Fp3(t[0],t[0]);
_Mul_Fp3(t[1],t[0],t[1]);
_Pow_Fp3(t[1],Value.Field.p_1_div_2,t[0]);
result:=t[0].IsOne;
end;
       {****** Test if an Fp3 is a Null *******}
function _Is_Fp3_Null(const value:Fp3Int):boolean;
begin
Result:=(value.a.Data.i16[-2]=0) and (Value.b.Data.i16[-2]=0)and (Value.c.Data.i16[-2]=0);
end;



{*******************************************************************************}
      ///  Definitions of an Fp3 integer operators and functions
{*******************************************************************************}
class operator Fp3Int.Add(const Left, Right: Fp3Int): Fp3Int;
begin
 if Left.Field.p<>Right.Field.p then raise Exception.Create('Ne peut additionner deux éléments de Fp3 avec des Modulo différents')
 else _Add_Fp3(Left,Right,Result);
end;

class operator  Fp3Int.Add(const Left:LInt; Right: Fp3Int): Fp3Int;
begin
_Add_FP_Fp3(Left,Right,Result);
end;

class operator  Fp3Int.Add(const Left: Fp3Int;Right:LInt): Fp3Int;
begin
_Add_Fp3_FP(Left,Right,Result);
end;

class operator Fp3Int.Equal(const Left: Fp3Int; Right: LInt): Boolean;
begin
Result:=_Compare_Fp3_FP(Left,Right);
end;

class operator Fp3Int.Equal(const Left, Right: Fp3Int): Boolean;
begin
Result:=_Compare_Fp3(Left,Right)=0;
end;

class operator Fp3Int.GreaterThan(const Left, Right: Fp3Int): Boolean;
begin
result:=_Compare_Fp3(Left,Right)=1;
end;

class operator Fp3Int.GreaterThanOrEqual(const Left,
  Right: Fp3Int): Boolean;
begin
Result:=_Compare_Fp3(Left,Right)>=0;
end;

class operator Fp3Int.LessThan(const Left, Right: Fp3Int): Boolean;
begin
Result:=_Compare_Fp3(Left,Right)=-1;
end;

class operator Fp3Int.LessThanOrEqual(const Left, Right: Fp3Int): Boolean;
begin
Result:=_Compare_Fp3(Left,Right)<=0;
end;

class operator Fp3Int.Subtract(const Left:LInt; Right: Fp3Int): Fp3Int;
begin
_Sub_FP_Fp3(Left,Right,Result);
end;

class operator Fp3Int.Subtract(const Left:Fp3Int; Right: LInt): Fp3Int;
begin
_Sub_Fp3_FP(Left,Right,Result);
end;

function Fp3Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.BitLength mod 8) * 8;
Setlength(Result,size*3);
Move(a.Data.i8[0],Result[0],size);
Move(b.Data.i8[0],Result[size],size);
Move(c.Data.i8[0],Result[2*size],size);
end;

class operator Fp3Int.Subtract(const Left, Right: Fp3Int): Fp3Int;
begin
_Sub_Fp3(Left,Right,Result);
end;

class operator Fp3Int.Multiply(const Left, Right: Fp3Int): Fp3Int;
begin
  if Left.Field.p<>Right.Field.p then raise Exception.Create('Ne peut multiplier deux éléments de Fp3 avec des Modulo différents')
  else _Mul_Fp3(Left,Right,Result);
end;

class operator Fp3Int.Multiply(const Left: Fp3Int; Right: Word): Fp3Int;
begin
_Mul_Fp3_Word(Left,Right,Result);
end;

class operator Fp3Int.Multiply(Left: Word; const Right: Fp3Int): Fp3Int;
begin
_Mul_Word_Fp3(Left,Right,Result);
end;

class operator Fp3Int.Multiply(Left: LInt; const Right: Fp3Int): Fp3Int;
begin
_Mul_FP_Fp3(Left,Right,Result);
end;

class operator Fp3Int.Negative(const Value: Fp3Int): Fp3Int;
begin
_Neg_Fp3(value,Result);
end;

class operator Fp3Int.NotEqual(const Left, Right: Fp3Int): Boolean;
begin
result:=_Compare_Fp3(Left,Right)<>0;
end;

class operator Fp3Int.Modulus(const Left: Fp3Int; Right: LInt): Fp3Int;
begin
_Mod_Fp3_FP(Left,Right,Result);
end;

function Fp3Int.Inverse: Fp3Int;
begin
_Inv_Fp3(Self,Result);
end;

procedure Fp3Int.SetFieldParams(FieldParams:PtrFieldParams);
begin
Field:=FieldParams;
end;

procedure Fp3Int.SetToRandom;
begin
Randomize;
GetRandomLIntLowerThan(a,Field^.p);
GetRandomLIntLowerThan(b,Field^.p);
GetRandomLIntLowerThan(c,Field^.p);
end;

procedure Fp3Int.SetToZero;
begin
a:=0;
b:=0;
c:=0;
end;

procedure Fp3Int.SetToOne;
begin
a:=1;
b:=0;
c:=0;
end;

function Fp3Int.Sqr: Fp3Int;
begin
_Sqr_Fp3(Self,Result);
end;

function Fp3Int.Sqrt:Fp3Int;
begin
_Sqrt_Fp3(Self,Result);
end;

function Fp3Int.Pow(Exponent : LInt): Fp3Int;
begin
_Pow_Fp3(Self,Exponent,Result);
end;

function Fp3Int.IsASquare: boolean;
begin
Result:=_Is_Fp3_ASquare(Self);
end;

function Fp3Int.IsOne: boolean;
begin
Result:=(a=1)and(b=0)and(c=0);
end;

function Fp3Int.IsZero: boolean;
begin
Result:=(a=0)and(b=0)and(c=0);
end;

function Fp3Int.toHexString: string;
begin
Result:=_Fp3_To_HexString(Self)
end;

function Fp3Int.toPowerP: Fp3Int;
begin
_PowP_Fp3(Self,Result);
end;

function Fp3Int.toDecimalString: string;
begin
Result:=_Fp3_To_DecimalString(Self)
end;

procedure Fp3Int.SetFormString(const Value: string);
begin
_Fp3_From_String(Value,Self);
end;

end.

