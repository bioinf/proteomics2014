#pragma once

#include <string>

struct Sequence {
    Sequence(std::string const &id, std::string const &sequence) :
        id(id),
        sequence(sequence) { }

    std::string id;
    std::string sequence;
};
