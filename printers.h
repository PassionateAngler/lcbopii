#ifndef _PRINTERS_H
#define _PRINTERS_H
#include <iostream>
#include <stack>
#include <sstream>
#include "atom.h"

using namespace simul;
class SdfPrinter:public AtomWalker
{
public:
	std::string get_sdf(Atom * start);
	/*void print();
	 SdfPrinter(Atom * start);
	*/

private:
	std::stack<Atom *> atoms;
	std::stack<std::pair<int,int>> bonds;
	static const std::string INDENT;

	virtual void bound_found(Atom *u, Atom *v)
	{
		bonds.push({(u->get_id())+1, (v->get_id())+1});
	};
	virtual void atom_found(Atom *u)
	{
		atoms.push(u);
	};

	std::string dump_atoms();
	std::string dump_bonds();
};

#endif
