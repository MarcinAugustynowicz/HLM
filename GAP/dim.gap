Read("./RatProbAlgTori/cryst3.gap");
Read("./RatProbAlgTori/cryst4.gap");
Read("./RatProbAlgTori/cryst5.gap");
Read("./RatProbAlgTori/cryst6.gap");
LoadPackage("Sonata");

crystQ := [cryst3, cryst4, cryst5, cryst6];

for dim in [4, 5, 6] do
	
	#lists all subgroups of GL(dim - 1, Z) and all subgroups of products of symmetric groups of the desired form.
	partitions := Partitions(dim);
	subgroups := [];
	for p in partitions do
		syms := List(p + 1, n -> SymmetricGroup(n));
		res := DirectProductOp(syms, syms[1]);
		Add(subgroups, res);
	od;
	Apply(subgroups, x -> Subgroups(x));
	subgroups2 := Concatenation(subgroups);
	Append(subgroups2, crystQ[dim-3]);
	cryst := crystQ[dim-2];
	i := 0;
	#for each finite subgroup of GL(dim - 1, Z), we check if it's isomorphic to some of the groups predicted by our conjecture. If it's not, then all Z-conjugacy classes of a group and lists of corresponding generators are written to files, for further processing by a Python script.
	for g in cryst do
		i := i + 1;
		iso := 0;
		for g2 in subgroups2 do
			if IsIsomorphicGroup(g, g2) then
				iso := 1;
				break;
			fi;
		od;
		if iso = 0 then
			crystZ := ZClassRepsQClass(g);
			for j in [1..Length(crystZ)] do
				name := StringFormatted("./dim{}/{}-{}.grp", dim, i, j);
				namegens := StringFormatted("./dim{}/{}-{}.gens", dim, i, j);
				PrintTo(name, "");
				PrintTo(namegens, "");
				for gg in crystZ[j] do					
					AppendTo(name, String(gg));
				od;	
				for gen in GeneratorsOfGroup(crystZ[j]) do
					AppendTo(namegens, String(gen));
				od;
			od;
		fi;
	od;
od;
