#include "Tree.h"

Tree::Tree(NodeNum num) {
    parent_ = nullptr;
    num_ = num;
}

Tree::Tree() {}

Tree::~Tree() {
    for (auto* child : children_) {
        delete child;
    }
}

Tree::const_iterator Tree::begin() { 
    return children_.begin(); 
}

Tree::const_iterator Tree::end() { 
    return children_.end(); 
}

void Tree::add_child(Tree * node) {
    children_.push_back(node);
    node->set_parent(this);
}

void Tree::remove_child(Tree * node) {
    list<Tree * >::iterator it = children_.begin();
    while (it != children_.end()) {
        if (*it == node) {
            children_.remove(*(it++));
        } else {
            ++it;
        }
    }
}

Tree * Tree::find(NodeNum node_number) {
    for (auto child : children_) {
        if (child->num_ == node_number) {
            return child;
        }
        child->find(node_number);
    }
    return nullptr;
}

Tree * Tree::get_parent() {
    return parent_;
}

void Tree::set_parent(Tree * parent) {
    parent_ = parent;
}

NodeNum Tree::get_num() {
    return num_;
}

void Tree::set_num(NodeNum num) {
    num_ = num;
}

void Tree::print_tree() {
    this->print(0);
}

void Tree::print(size_t space_count) {
    if (this != nullptr) {
        cout << string(space_count, ' ') << this->num_ << endl;
        space_count += 2;
        for (auto* child : this->children_) {
            child->print(space_count);
        }
    }
}

