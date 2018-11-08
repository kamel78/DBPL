unit Fp2Arithmetic;

interface
uses System.SysUtils,VCL.dialogs, LargeIntegers;

type

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp2. Fp2 is the quadratic extention of Fp with respect to the
   irreducible polynômial u²+Beta (u^2=Beta). Elements are on the form a+u.b (a,b from Fp).

********************************************************************************************}

             { Definition of the FP2 field Extension of the field FP modulo the irriductibe polynomial X^2+Beta=0 }

  Fp2Int=record
  Field:PtrFieldParams;
  a,b:LInt;
    public
    function toHexString: string;
    function toDecimalString: string;
    function toByteArray:Tbytes;
    function Inverse:Fp2Int;
    function Sqr:Fp2Int;
    function Sqrt:Fp2Int;
    function Pow(Exponent: LInt): Fp2Int;
    function Conjugate:Fp2Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function IsASquare:boolean;
    procedure SetToRandom;
    procedure SetToZero;
    procedure SetToOne;
    procedure SetFieldParams(FieldParams:PtrFieldParams);
    procedure SetFormString(const Value: string);
    class operator Add(const Left, Right: Fp2Int): Fp2Int;
    class operator Add(const Left:LInt; Right: Fp2Int): Fp2Int;
    class operator Add(const Left: Fp2Int;Right:LInt): Fp2Int;
    class operator Subtract(const Left, Right: Fp2Int): Fp2Int;
    class operator Subtract(const Left:LInt; Right: Fp2Int): Fp2Int;
    class operator Subtract(const Left:Fp2Int; Right: LInt): Fp2Int;
    class operator Multiply(const Left, Right: Fp2Int): Fp2Int;
    class operator Multiply(const Left: Fp2Int; Right: Word): Fp2Int;
    class operator Multiply(Left: Word; const Right: Fp2Int): Fp2Int;
    class operator Multiply(Left: LInt; const Right: Fp2Int): Fp2Int;
    class operator Modulus(const Left: Fp2Int; Right: LInt): Fp2Int;
    class operator Negative(const Value: Fp2Int): Fp2Int;
    class operator Equal(const Left, Right: Fp2Int): Boolean;
    class operator Equal(const Left:Fp2Int;Right:LInt): Boolean;
    class operator NotEqual(const Left, Right: Fp2Int): Boolean;
    class operator GreaterThan(const Left, Right: Fp2Int): Boolean;
    class operator GreaterThanOrEqual(const Left, Right: Fp2Int): Boolean;
    class operator LessThan(const Left, Right: Fp2Int): Boolean;
    class operator LessThanOrEqual(const Left, Right: Fp2Int): Boolean;
  end;

  procedure _Add_FP2(Left,Right:Fp2Int; var Result:Fp2Int);
  procedure _Add_FP_FP2(Left:LInt; Right: Fp2Int; var Result:Fp2Int);
  procedure _Add_FP2_FP(Left:Fp2Int; Right: LInt; var Result:Fp2Int);
  procedure _Sub_FP2(Left,Right:Fp2Int;var Result: Fp2Int);
  procedure _Sub_FP2_FP(Left:Fp2Int; Right: LInt; var Result:Fp2Int);
  procedure _Conjugate_FP2(const Source:Fp2Int;var result:Fp2Int);
  function  _Compare_FP2(const Left,Right:Fp2Int):shortint;
  function  _Compare_FP2_FP(const Left:Fp2Int;Right:LInt):boolean;
  procedure _Mul_FP2(const Left,Right:Fp2Int;var Result:Fp2Int);
  procedure _Mul_FP2_FP(const Left: Fp2Int; Right: LInt; var Result: Fp2Int);
  procedure _Mul_FP_FP2(const Left: LInt; Right: Fp2Int; var Result: Fp2Int);
  procedure _Mul_FP2_Word(const Left: Fp2Int; Right: Word; var Result: Fp2Int);
  procedure _Mul_Word_FP2(const Left: Word; Right: Fp2Int; var Result: Fp2Int);
  procedure _Neg_FP2(const Value: Fp2Int; var Result: Fp2Int);
  procedure _Mod_FP2_FP(const Left: Fp2Int; Right: LInt; var Result: Fp2Int);
  procedure _Sqr_FP2(const Value:Fp2Int;var Result:Fp2Int);
  procedure _Inv_FP2(const value:Fp2Int;var Result:Fp2Int);
  procedure _Sqrt_FP2(const Value:Fp2Int;var Result:Fp2Int);
  procedure _Pow_FP2(const value:Fp2Int;Exponent:LInt;var Result:Fp2Int);
  function  _FP2_To_HexString(const value:Fp2Int):string;
  procedure _FP2_From_String(const value:String;var Result:Fp2Int);
  function  _Is_FP2_ASquare(const value:Fp2Int):boolean;
  function  _Is_FP2_Null(const value:Fp2Int):boolean;
  function  InitFieldParams(Beta:integer;P:Lint):PtrFieldParams;
  procedure _NoMod_Add_FP2(Left,Right:Fp2Int; var Result:Fp2Int);
  procedure _NoMod_Sub_FP2(Left,Right:Fp2Int;var Result: Fp2Int);
  procedure _Nomod_Sqr_FP2(const Value:Fp2Int;var Result:Fp2Int);
  procedure _Nomod_Mul_FP2(const Left,Right:Fp2Int;var Result:Fp2Int);

implementation

function InitFieldParams(Beta:integer;P:Lint):PtrFieldParams;
begin
New(Result);
Result.Beta:=Beta;
Result^.p:=p;
Result^.p_1_div_2:=(p-1) shr 1;
Result^.inv_2:=LInt(2).InversModulo(p);
New(Result.MontgomeryData);
Result.pmod8:=p mod 8;
InitMontgomeryStruct(P,Result.MontgomeryData);

/// here is for compatibility in Fp6 for KSS36 (Sqrt especially)
Result^.sru:=Lint(Beta).PowMod((p-1)/6,p,Result.MontgomeryData); // For Square Root Computation
_Sqr_LInt(Result.sru,Result.sru2);
_Mod_LInt(Result.sru2,p,Result.sru2);
Result.pmod8:=p mod 8;
Result.p_5_div_8:=(p-5)/8;
Result.p_3_div_4:=(p-3)/4;

end;

{*******************************************************************************}
            ///      Procedures for Fp2 Arithmetic
{*******************************************************************************}

        {**********   Add two FP2 integers *****************}
        // Without Modular reduction (for optimization in some cases)
procedure _NoMod_Add_FP2(Left,Right:Fp2Int; var Result:Fp2Int);
begin
Result.Field:=Left.Field;
_Add_LInt(Left.a,Right.a,Result.a);
_Add_LInt(Left.b,Right.b,Result.b);
end;

        {**********   Add two FP2 integers *****************}
procedure _Add_FP2(Left,Right:Fp2Int; var Result:Fp2Int);
begin
Result.Field:=Left.Field;
_Add_LInt(Left.a,Right.a,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Add_LInt(Left.b,Right.b,Result.b);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
end;
        {**********   Add LInt with FP2 integer***********}
procedure _Add_FP_FP2(Left:LInt; Right: Fp2Int; var Result:Fp2Int);
begin
Result.Field:=Right.Field;
_Add_LInt(Right.a,Left,Result.a);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
Result.b:=Right.b;
end;
        {**********   Add FP2 integer with LInt **********}
procedure _Add_FP2_FP(Left:Fp2Int; Right: LInt; var Result:Fp2Int);
begin
Result.Field:=Left.Field;
_Add_LInt(Left.a,Right,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
Result.b:=Left.b;
end;
        {**********   Substract two FP2 integers ***********}
        // Without Modular reduction (for optimization in some cases)
procedure _NoMod_Sub_FP2(Left,Right:Fp2Int;var Result: Fp2Int);
begin
Result.Field:=Left.Field;
_Sub_LInt(Left.a,Right.a,Result.a);
_Sub_LInt(Left.b,Right.b,Result.b);
end;

        {**********   Substract two FP2 integers ***********}
procedure _Sub_FP2(Left,Right:Fp2Int;var Result: Fp2Int);
begin
Result.Field:=Left.Field;
_Sub_LInt(Left.a,Right.a,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Sub_LInt(Left.b,Right.b,Result.b);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
end;
        {*********   Substract LInt with FP2 integer *****}
procedure _Sub_FP_FP2(const Left:LInt; Right: Fp2Int; var Result:Fp2Int);
begin
Result.Field:=Right.Field;
Result.b:=Right.b;
_Sub_LInt(Left,Right.a,Result.a);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
end;
        {******** Substract FP2 integer with LInt ********}
procedure _Sub_FP2_FP(Left:Fp2Int; Right: LInt; var Result:Fp2Int);
begin
Result.Field:=Left.Field;
Result.b:=Left.b;
_Sub_LInt(Left.a,Right,Result.a);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
end;
        {**********   Conjugate an FP2 integer *************}
procedure _Conjugate_FP2(const Source:Fp2Int;var result:Fp2Int);
begin
Result.Field:=Source.Field;
Result.a:=Source.a;
if Source.b<>0 then _Sub_LInt(Source.Field^.p,Source.b,Result.b) else Result.b:=0;
end;
        {**********   Compare two FP2 integers *************}
function _Compare_FP2(const Left,Right:Fp2Int):shortint;
begin
if (_Compare_LInt(Left.a,right.a)=0) and(_Compare_LInt(left.b,right.b)=0) then result:=0
else result:=-1;
end;
        {******** Compare FP2 integer with LInt **********}
function _Compare_FP2_FP(const Left:Fp2Int;Right:LInt):boolean;
begin
Result:=(Left.b.Data.i16[-2]=0)and(_Compare_LInt(Left.a,Right)=0);
end;

        {********** Multiply two FP2 integers *****************}
        // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Mul_FP2(const Left,Right:Fp2Int;var Result:Fp2Int);
var t:array[0..3] of LInt;
begin
Result.Field:=Left.Field;
if (Right.a.Data.i16[-2]=1)
   and(Right.b.Data.i16[-2]=1)
   and(Right.a.Data.i16[0]=1)
   and(Right.b.Data.i16[0]=1)then begin
                                  // Optimization if sigma=1+u (generally the case)
                                  _Sub_Lint(Left.a,Left.b,Result.a);
                                  _Add_Lint(Left.a,Left.b,Result.b);
                                  end
else begin
      // with Karatsuba Multiplication
     _Add_LInt(Left.a,Left.b,t[3]);
     _Add_LInt(Right.a,Right.b,t[0]);
     _Mul_LInt(t[3],t[0],t[1]);
     _Mul_LInt(Left.a,Right.a,t[2]);
     _Sub_LInt(t[1],t[2],t[3]);
     _Mul_LInt(Left.b,Right.b,t[1]);
     _Sub_LInt(t[3],t[1],Result.b);
     if Left.Field^.Beta=-1 then begin
                                 t[1].Data.i16[-1]:=not t[1].Data.i16[-1];
                                 _Add_LInt(t[2],t[1],Result.a);
                                 end
     else begin
          _Mul_LInt(t[1],LInt(Left.Field^.Beta),t[0]);
          _Add_LInt(t[2],t[0],Result.a);
          end;
     end;
end;


        {********** Multiply two FP2 integers *****************}
procedure _Mul_FP2(const Left,Right:Fp2Int;var Result:Fp2Int);
var t:array[0..3] of LInt;
begin
Result.Field:=Left.Field;
if (Right.a.Data.i16[-2]=1)and(Left.Field.Beta=-1)
   and(Right.b.Data.i16[-2]=1)
   and(Right.a.Data.i16[0]=1)
   and(Right.b.Data.i16[0]=1)then begin
                                  // Optimization if sigma=1+u (generally the case)
                                  _Sub_Lint(Left.a,Left.b,t[2]);
                                  _Add_Lint(Left.a,Left.b,t[3]);
                                  end
else begin
      // with Karatsuba Multiplication
     _Add_LInt(Left.a,Left.b,t[3]);
     _Add_LInt(Right.a,Right.b,t[0]);
     _Mul_LInt(t[3],t[0],t[1]);
     _Mul_LInt(Left.a,Right.a,t[2]);
     _Sub_LInt(t[1],t[2],t[3]);
     _Mul_LInt(Left.b,Right.b,t[1]);
     _Sub_LInt(t[3],t[1],t[3]);
     if Left.Field^.Beta=-1 then begin
                                 t[1].Data.i16[-1]:=not t[1].Data.i16[-1];
                                 _Add_LInt(t[2],t[1],t[2]);
                                 end
     else begin
          _Mul_LInt(t[1],LInt(Left.Field^.Beta),t[0]);
          _Add_LInt(t[2],t[0],t[2]);
          end;
     end;
_Mod_LInt(t[2],Left.Field^.p,Result.a);
_Mod_LInt(t[3],Left.Field^.p,Result.b);
end;
        {********** Multiply FP2 with LInt *****************}
procedure _Mul_FP2_FP(const Left: Fp2Int; Right: LInt; var Result: Fp2Int);
begin
Result.Field:=Left.Field;
_Mul_LInt(Left.a,Right,Result.a);
_Mul_LInt(Left.b,Right,Result.b);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
end;
        {********** Multiply LInt with FP2 *****************}
procedure _Mul_FP_FP2(const Left: LInt; Right: Fp2Int; var Result: Fp2Int);
begin
Result.Field:=Right.Field;
_Mul_LInt(Left,Right.a,Result.a);
_Mul_LInt(Left,Right.b,Result.b);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
_Mod_LInt(Result.b,Right.Field^.p,Result.b);
end;
        {********** Multiply FP2 with Word *****************}
procedure _Mul_FP2_Word(const Left: Fp2Int; Right: Word; var Result: Fp2Int);
begin
Result.Field:=Left.Field;
_Mul_LInt(Left.a,LInt(Right),Result.a);
_Mul_LInt(Left.b,LInt(Right),Result.b);
_Mod_LInt(Result.a,Left.Field^.p,Result.a);
_Mod_LInt(Result.b,Left.Field^.p,Result.b);
end;
        {********** Multiply Word with FP2 *****************}
procedure _Mul_Word_FP2(const Left: Word; Right: Fp2Int; var Result: Fp2Int);
begin
Result.Field:=Right.Field;
_Mul_LInt(LInt(Left),Right.a,Result.a);
_Mul_LInt(LInt(Left),Right.b,Result.b);
_Mod_LInt(Result.a,Right.Field^.p,Result.a);
_Mod_LInt(Result.b,Right.Field^.p,Result.b);
end;
        {********** Get negative of an FP2 *****************}
procedure _Neg_FP2(const Value: Fp2Int; var Result: Fp2Int);
begin
result.Field:=Value.Field;
if Value.a.Data.i16[-2]<>0 then  _Sub_LInt(Value.Field^.p,Value.a,Result.a) else Result.a.Data.i16[-2]:=0;
if Value.b.Data.i16[-2]<>0 then  _Sub_LInt(Value.Field^.p,Value.b,Result.b)else Result.b.Data.i16[-2]:=0;
end;
        {********** Get modulo of an FP2 *****************}
procedure _Mod_FP2_FP(const Left: Fp2Int; Right: LInt; var Result: Fp2Int);
begin
Result.Field:=Left.Field;
_Mod_LInt(Left.a,Right,Result.a);
_Mod_LInt(Left.b,Right,Result.b);
end;

        {********** Get Square of an FP2 *****************}
        // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Sqr_FP2(const Value:Fp2Int;var Result:Fp2Int);
var t:array[0..2]of LInt;
begin
// with Karatsuba Squaring
Result.Field:=Value.Field;
if Value.Field.Beta=-1 then begin
                            _Add_Lint(value.a,value.b,t[0]);
                            _Sub_LInt(value.a,value.b,t[1]);
                            _Mul_LInt(t[0],t[1],t[2]);
                            _Add_Lint(value.b,value.b,t[0]);
                            _Mul_LInt(t[0],value.a,Result.b);
                            _HCopy_LInt(t[2],Result.a);
                            end
else begin
     _Sqr_LInt(Value.a,t[0]);
     _Sqr_LInt(Value.b,t[1]);
     _Add_LInt(Value.a,Value.b,Result.b);
     _Sqr_LInt(Result.b,Result.a);
     _Sub_LInt(Result.a,t[0],Result.b);
     _Sub_LInt(Result.b,t[1],Result.b);
     _Mul_LInt(t[1],LInt(Value.Field^.beta),Result.a);
     _Add_LInt(Result.a,t[0],Result.a);
     end;
end;

        {********** Get Square of an FP2 *****************}
procedure _Sqr_FP2(const Value:Fp2Int;var Result:Fp2Int);
var t:array[0..2]of LInt;
begin
// with Karatsuba Squaring
Result.Field:=Value.Field;
if Value.Field.Beta=-1 then begin
                            _Add_Lint(value.a,value.b,t[0]);
                            _Sub_LInt(value.a,value.b,t[1]);
                            _Mul_LInt(t[0],t[1],t[2]);
                            _Add_Lint(value.b,value.b,t[0]);
                            _Mul_LInt(t[0],value.a,Result.b);
                            _Mod_LInt(t[2],Value.Field^.p,Result.a);
                            _Mod_LInt(Result.b,Value.Field^.p,Result.b);
                            end
else begin
     _Sqr_LInt(Value.a,t[0]);
     _Sqr_LInt(Value.b,t[1]);
     _Add_LInt(Value.a,Value.b,Result.b);
     _Sqr_LInt(Result.b,Result.a);
     _Sub_LInt(Result.a,t[0],Result.b);
     _Sub_LInt(Result.b,t[1],Result.b);
     _Mul_LInt(t[1],LInt(Value.Field^.beta),Result.a);
     _Add_LInt(Result.a,t[0],Result.a);
     _Mod_LInt(Result.a,Value.Field^.p,Result.a);
     _Mod_LInt(Result.b,Value.Field^.p,Result.b);
     end;
end;
        {********** Get inverse of an FP2 *****************}
procedure _Inv_FP2(const value:Fp2Int;var Result:Fp2Int);
var t:array[0..2] of LInt;
begin
Result.Field:=Value.Field;
_Sqr_LInt(Value.b,t[0]);
_Mul_LInt(t[0],LInt(Value.Field^.beta),t[1]);
_Sqr_LInt(Value.a,t[0]);
_Sub_LInt(t[0],t[1],t[0]);
_Inv_Mod_LInt(t[0],Value.Field^.p,t[1]);
_Mul_LInt(Value.a,t[1],t[2]);
_Sub_LInt(Value.Field^.p,Value.b,t[0]);
_Mul_LInt(t[0],t[1],Result.b);
_Mod_LInt(t[2],Value.Field^.p,Result.a);
_Mod_LInt(Result.b,Value.Field^.p,Result.b);
end;
        {********** Get Square root of an FP2 ************}
procedure _Sqrt_FP2(const Value:Fp2Int;var Result:Fp2Int);
var t:array[0..2] of LInt;
    delta:LInt;
begin
Result.Field:=Value.Field;
_Sqr_LInt(Value.b,t[0]);
_Mul_LInt(t[0],LInt(Value.Field^.beta),t[1]);
_Sqr_LInt(Value.a,Delta);
_Sub_LInt(Delta,t[1],Delta);
_Mod_LInt(Delta,Value.Field^.p,Delta);
if not ModSquareLInt(delta,Value.Field^.p,t[0],Value.Field.MontgomeryData, False) then
  raise ERangeError.Create( 'This FP2 element is not a Square .....')
else begin
     _Add_LInt(Value.a,t[0],t[2]);
     _Mul_LInt(t[2],Value.Field^.inv_2,t[1]);
     _Mod_LInt(t[1],Value.Field^.p,t[1]);
     if not ModSquareLInt(t[1],Value.Field^.p,Result.a,Value.Field.MontgomeryData,False) then
                  begin
                  _Sub_LInt(t[1],t[0],t[1]);
                  _Mod_LInt(t[1],Value.Field^.p,t[1]);
                  if not ModSquareLInt(t[1],Value.Field^.p,Result.a,Value.Field.MontgomeryData,False) then
                          begin
                          _Sub_LInt(Value.Field^.p,t[1],t[1]);
                          _Mod_LInt(t[1],Value.Field^.p,t[1]);
                          if not ModSquareLInt(t[1],Value.Field^.p,Result.a,Value.Field.MontgomeryData,False) then
                                  begin
                                  _Sub_LInt(t[1],t[0],t[1]);
                                  _Mod_LInt(t[1],Value.Field^.p,t[1]);
                                  ModSquareLInt(t[1],Value.Field^.p,Result.a,Value.Field.MontgomeryData,False)
                                  end;
                          end;
                  end;
     if (Result.a.Data.i16[-2]=0) then Result.b.Data.i32[-1]:=0
     else begin
          _Hcopy_LInt(Result.a,t[0]);
          _Shl_LInt(t[0],1);
          _Mod_LInt(t[0],Value.Field^.p,t[0]);
          _Inv_Mod_LInt(t[0],Value.Field^.p,t[1]);
          _Mul_LInt(Value.b,t[1],Result.b);
          _Mod_LInt(Result.b,Value.Field^.p,Result.b);
          end;
     end;
end;
        {********** Raise an FP2 to a LInt power ************}
procedure _Pow_FP2(const value:Fp2Int;Exponent:LInt;var Result:Fp2Int);
var i:word;
begin
Result.Field:= Value.Field;
if (Value.b.Data.i16[-2]=0) then begin
                                 _Mg_Mod_Pow_LInt(Value.a,Exponent,Value.Field^.p,Value.Field.MontgomeryData,Result.a);
                                 Result.b.Data.i16[-2]:=0;
                                 end
else begin
     Result.a:=Value.a;
     Result.b:=Value.b;
     if (Exponent>1) then
     for i:=Exponent.BitLength-2 downto 0 do begin
                                              _Sqr_FP2(Result,Result);
                                              if _Is_BitSet_At(Exponent,i) then
                                              _Mul_FP2(Result,Value,Result);
                                              end;
     end;
end;

        {****** Convert an FP2 to a Decimal String **********}
function _FP2_To_DecimalString(const value:Fp2Int):string;
begin
  if (Value.a.Data.i16[-2]=0)then begin
                   if (Value.b.Data.i16[-2]=0)then Result:='0'
                   else begin
                        if Value.b=1 then Result:='u'
                        else Result:='u * '+Value.b.ToDecimalString;
                        end;
                   end
  else begin
       if (Value.b.Data.i16[-2]=0) then Result:=Value.a.ToDecimalString
       else begin
            if Value.b=1 then Result:='u +'+Value.a.ToDecimalString
            else Result:='u * '+Value.b.ToDecimalString+' + '+Value.a.ToDecimalString;
            end;
       end;
end;

        {****** Convert an FP2 to a Hexadecimal String *******}
function _FP2_To_HexString(const value:Fp2Int):string;
begin
  if (Value.a.Data.i16[-2]=0)then begin
                   if (Value.b.Data.i16[-2]=0) then Result:='0'
                   else begin
                        if Value.b=1 then Result:='u'
                        else Result:='u * '+Value.b.ToHexString;
                        end;
                   end
  else begin
       if (Value.b.Data.i16[-2]=0) then Result:=Value.a.ToHexString
       else begin
            if Value.b=1 then Result:='u +'+Value.a.ToHexString
            else Result:='u * '+Value.b.ToHexString+' + '+Value.a.ToHexString;
            end;
       end;
end;
        {****** Convert String to an FP2 (Decimal/Hex)*******}
procedure _FP2_From_String(const value:String;var Result:Fp2Int);
var i:integer;
    s1,s2,Val:string;
    valid:boolean;
begin
Val:=value;
if (pos('+u*',Val)=0) then begin
                            Valid:=true;
                            for i:=1 to length(Val) do if not(Val[i]in ['x','-','$','0'..'9','a'..'f','A'..'F']) then valid:=false;
                            if not valid then raise Exception.Create('Valeur Invalide pout un type FP2Int')
                            else begin
                                 Result.a:=Val;
                                 Result.b:=0;
                                 end;
                            end
else begin
     s1:=copy(Val,1,pos('+u*',Val)-1);
     s2:=copy(Val,pos('+u*',Val)+3,length(Val));
     Valid:=true;
     for i:=1 to length(s1) do if not(S1[i]in ['x','-','$','0'..'9','a'..'f','A'..'F']) then valid:=false;
     if Valid then begin
                   for i:=1 to length(s2) do if not(S2[i]in ['x','-','$','0'..'9','a'..'f','A'..'F']) then valid:=false;
                   if Valid then Begin
                                 Result.a:=S1;
                                 Result.b:=s2;
                                 end
                   else raise Exception.Create('Valeur Invalide pout un type FP2Int');
                   end;
     end;
end;

       {****** Test if an FP2 is a Sequare *******}
function _Is_FP2_ASquare(const value:Fp2Int):boolean;
var t:array[0..1] of LInt;
begin
_Sqr_LInt(Value.b,t[0]);
_Mul_LInt(LInt(Value.Field^.beta),t[0],t[1]);
_Sqr_LInt(Value.a,t[0]);
_Sub_LInt(t[0],t[1],t[0]);
_Mod_LInt(t[0],Value.Field^.p,t[0]);
_Mg_Mod_Pow_LInt(t[0],Value.Field^.p_1_div_2,Value.Field^.p,Value.Field.MontgomeryData,t[1]);
Result:=(_Compare_LInt(t[1],1)=0);
end;
       {****** Test if an FP2 is a Null *******}
function _Is_FP2_Null(const value:Fp2Int):boolean;
begin
Result:=(value.a.Data.i16[-2]=0) and (Value.b.Data.i16[-2]=0);
end;



{*******************************************************************************}
      ///  Definitions of an FP2 integer operators and functions
{*******************************************************************************}
class operator Fp2Int.Add(const Left, Right: Fp2Int): Fp2Int;
begin
 if Left.Field.p<>Right.Field.p then raise Exception.Create('Ne peut additionner deux éléments de FP2 avec des Modulo différents')
 else _Add_FP2(Left,Right,Result);
end;

class operator  Fp2Int.Add(const Left:LInt; Right: Fp2Int): Fp2Int;
begin
_Add_FP_FP2(Left,Right,Result);
end;

class operator  Fp2Int.Add(const Left: Fp2Int;Right:LInt): Fp2Int;
begin
_Add_FP2_FP(Left,Right,Result);
end;

function Fp2Int.Conjugate: Fp2Int;
begin
_Conjugate_FP2(Self,Result);
end;

class operator Fp2Int.Equal(const Left: Fp2Int; Right: LInt): Boolean;
begin
Result:=_Compare_FP2_FP(Left,Right);
end;

class operator Fp2Int.Equal(const Left, Right: Fp2Int): Boolean;
begin
Result:=_Compare_FP2(Left,Right)=0;
end;

class operator Fp2Int.GreaterThan(const Left, Right: Fp2Int): Boolean;
begin
result:=_Compare_FP2(Left,Right)=1;
end;

class operator Fp2Int.GreaterThanOrEqual(const Left,
  Right: Fp2Int): Boolean;
begin
Result:=_Compare_FP2(Left,Right)>=0;
end;

class operator Fp2Int.LessThan(const Left, Right: Fp2Int): Boolean;
begin
Result:=_Compare_FP2(Left,Right)=-1;
end;

class operator Fp2Int.LessThanOrEqual(const Left, Right: Fp2Int): Boolean;
begin
Result:=_Compare_FP2(Left,Right)<=0;
end;

class operator Fp2Int.Subtract(const Left:LInt; Right: Fp2Int): Fp2Int;
begin
_Sub_FP_FP2(Left,Right,Result);
end;

class operator Fp2Int.Subtract(const Left:Fp2Int; Right: LInt): Fp2Int;
begin
_Sub_FP2_FP(Left,Right,Result);
end;

class operator Fp2Int.Subtract(const Left, Right: Fp2Int): Fp2Int;
begin
_Sub_FP2(Left,Right,Result);
end;

class operator Fp2Int.Multiply(const Left, Right: Fp2Int): Fp2Int;
begin
  if Left.Field.p<>Right.Field.p then raise Exception.Create('Ne peut multiplier deux éléments de FP2 avec des Modulo différents')
  else _Mul_FP2(Left,Right,Result);
end;

class operator Fp2Int.Multiply(const Left: Fp2Int; Right: Word): Fp2Int;
begin
_Mul_FP2_Word(Left,Right,Result);
end;

class operator Fp2Int.Multiply(Left: Word; const Right: Fp2Int): Fp2Int;
begin
_Mul_Word_FP2(Left,Right,Result);
end;

class operator Fp2Int.Multiply(Left: LInt; const Right: Fp2Int): Fp2Int;
begin
_Mul_FP_FP2(Left,Right,Result);
end;

class operator Fp2Int.Negative(const Value: Fp2Int): Fp2Int;
begin
_Neg_FP2(value,Result);
end;

class operator Fp2Int.NotEqual(const Left, Right: Fp2Int): Boolean;
begin
result:=_Compare_FP2(Left,Right)<>0;
end;

class operator Fp2Int.Modulus(const Left: Fp2Int; Right: LInt): Fp2Int;
begin
_Mod_FP2_FP(Left,Right,Result);
end;

function Fp2Int.Inverse: Fp2Int;
begin
_Inv_FP2(Self,Result);
end;

procedure Fp2Int.SetFieldParams(FieldParams:PtrFieldParams);
begin
Field:=FieldParams;
end;

procedure Fp2Int.SetToRandom;
begin
Randomize;
GetRandomLIntLowerThan(a,Field^.p);
GetRandomLIntLowerThan(b,Field^.p);
end;

procedure Fp2Int.SetToZero;
begin
a:=0;
b:=0;
end;

procedure Fp2Int.SetToOne;
begin
a:=1;
b:=0;
end;

function Fp2Int.Sqr: Fp2Int;
begin
_Sqr_FP2(Self,Result);
end;

function Fp2Int.Sqrt:Fp2Int;
begin
_Sqrt_FP2(Self,Result);
end;

function Fp2Int.Pow(Exponent : LInt): Fp2Int;
begin
_Pow_FP2(Self,Exponent,Result);
end;

function Fp2Int.IsASquare: boolean;
begin
Result:=_Is_FP2_ASquare(Self);
end;

function Fp2Int.IsOne: boolean;
begin
Result:=(a=1)and(b=0);
end;

function Fp2Int.IsZero: boolean;
begin
Result:=(a=0)and(b=0);
end;

function Fp2Int.toHexString: string;
begin
Result:=_FP2_To_HexString(Self)
end;

function Fp2Int.toByteArray: Tbytes;
var size:integer;
begin
size:=(a.BitLength mod 8) * 8;
Setlength(Result,size*2);
Move(a.Data.i8[0],Result[0],size);
Move(b.Data.i8[0],Result[size],size);
end;

function Fp2Int.toDecimalString: string;
begin
Result:=_FP2_To_DecimalString(Self)
end;

procedure Fp2Int.SetFormString(const Value: string);
begin
_FP2_From_String(Value,Self);
end;

end.
