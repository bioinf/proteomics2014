#include <iostream>
#include <list>
#include <string>

using std::cout;
using std::endl;
using std::list;
using std::string;

typedef int NodeNum;

struct Tree {
    typedef list<Tree *>::const_iterator const_iterator;
    Tree(NodeNum n);
    Tree();
    ~Tree();
    const_iterator begin();
    const_iterator end();
    void add_child(Tree * node);
    void remove_child(Tree * node);
    Tree * find(NodeNum node_number);
    Tree * get_parent();
    void set_parent(Tree * parent);
    NodeNum get_num();
    void set_num(NodeNum num);
    void print_tree();
    void print(size_t space_count);

    private:
        list<Tree *> children_;
        Tree * parent_;
        NodeNum num_;
};

