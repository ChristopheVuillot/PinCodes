CoxeterMatrix := function(S)
    #returns Coxeter matrix from Schlaefli symbol
    #works for general schlaefli symbol
    local C,D,i,j;
    D := Size(S);
    C := IdentityMat(D+1);
    for i in [1..D] do
        C[i][i+1] := S[i];
        C[i+1][i] := S[i];
    od;
    for i in [1..D+1] do
        for j in [i+2..D+1] do
            C[i][j] := 2;
            C[j][i] := 2;
        od;
    od;
    return C;
end;


CoxeterGroupFromMatrix := function(C)
  #given the Coxeter matrix C, returns the corresponding Coxeter group
  local D, F, gensF, rel, G, i, j;
  D := Size(C);
  F := FreeGroup(D);
  gensF := GeneratorsOfGroup(F);
  rel := [];
  for i in [1..D] do
	for j in [1..D] do
	    Add(rel,(gensF[i]*gensF[j])^C[i][j]);
	od;
    od;
    G:= F/rel;
    return G;
end;


Petri_compact_fromMatrix := function(C,r)
    #Given the Coxeter matrix C and a number r, returns the incidence polytope of type S
    #with Petri polygon - translation t^r added to the relations
    local G,Ggens,i,t,D,F,gensF,j,rel;
    D := Size(C);
    F := FreeGroup(D);
    gensF := GeneratorsOfGroup(F);

    rel := [];
    for i in [1..D] do
	for j in [1..D] do
	    Add(rel,(gensF[i]*gensF[j])^C[i][j]);
	od;
    od;

    t := Identity(F);
    for i in [1..Size(gensF)] do
	t := t*gensF[i];
    od;
    Add(rel,t^r);

    G:= F/rel;
    return G;
end;


Create_and_write_to_file_code := function(G, ggen, x, z, name)
  #Create a quantum codes from Coxeter group G with the group elements as qubits
  #the (n-x)-generators cosets as x-checks and (n-z)-generators cosets as z-checks
  #and write the two matrices in sparse format to a file.
  local n, xchecks, zchecks, numqubits, groupelements, numxchecks, numzchecks, filex, outputx, filez, outputz, eindex, lg, j, e, xmatrixlist, zmatrixlist;
  n := Length(ggen);
  xchecks := [];
  zchecks := [];
  for lg in Combinations(ggen, n-x) do
    Append(xchecks,RightCosets(G,Subgroup(G,lg)));
  od;
  for lg in Combinations(ggen, n-z) do
    Append(zchecks,RightCosets(G,Subgroup(G,lg)));
  od;
  numqubits := Size(G);
  groupelements := Elements(G);
  numxchecks := Size(xchecks);
  numzchecks := Size(zchecks);
  xmatrixlist := [];
  zmatrixlist := [];
  
  filex := Concatenation(["../PCMatrices/", name, "X.sms"]);
  outputx := OutputTextFile(filex, false );;
  AppendTo(outputx, Concatenation([String(numxchecks)," ",String(numqubits)," 2\n"]));
  for j in [1..numxchecks] do
    for e in Elements(xchecks[j]) do
      eindex := Position(groupelements, e);
      Append(xmatrixlist, [[j, eindex, 1]]);
      AppendTo(outputx, Concatenation([String(j)," ",String(eindex)," 1\n"]));
    od;
  od;
  AppendTo(outputx, "0 0 0\n");
  CloseStream(outputx);

  filez := Concatenation(["../PCMatrices/", name, "Z.sms"]);
  outputz := OutputTextFile(filez, false );;
  AppendTo(outputz, Concatenation([String(numzchecks)," ",String(numqubits)," 2\n"]));
  for j in [1..numzchecks] do
    for e in Elements(zchecks[j]) do
      eindex := Position(groupelements, e);
      Append(zmatrixlist, [[j, eindex, 1]]);
      AppendTo(outputz, Concatenation([String(j)," ",String(eindex)," 1\n"]));
    od;
  od;
  AppendTo(outputz, "0 0 0\n");
  CloseStream(outputz);

  return [numqubits, numxchecks, numzchecks, xmatrixlist, zmatrixlist];
end;


script := function(SchlaefliS, r, x, z, name)
  local G, ggen;
  G := Petri_compact_fromMatrix(CoxeterMatrix(SchlaefliS), r);
  ggen := GeneratorsOfGroup(G);
  Create_and_write_to_file_code(G, ggen, x, z, name);
end;

script([4,4,4,4], 2, 1, 2, "coxeter44442");
