LoadPackage("liepring", false);

#
PrintConditions := function(Ls)
  local i, c;
  for i in [1..Length(Ls)] do
    c:=LibraryConditions(Ls[i]);
    if c<>["",""] then
      Print(i, ": ", c, "\n");
    fi;
  od;
end;;

#
CountFamilies := function(Ls)
  local i, p, s, old;
  old:=InfoLevel(InfoPerformance);
  SetInfoLevel(InfoPerformance, 0);
  for i in [1..Length(Ls)] do
    Print(i, ": ");
    for p in [2,3,5,7,11,13] do
      s := LiePRingsInFamily(Ls[i], p);
      if s = fail then s := 0; else s := Length(s); fi;
      Print(s, ", ");
    od;
    Print("\n");
  od;
  SetInfoLevel(InfoPerformance, old);
end;;
