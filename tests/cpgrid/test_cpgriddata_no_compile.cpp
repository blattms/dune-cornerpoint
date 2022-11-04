#include "config.h"
#include<tuple>
#include<vector>
#include<memory>
#include <opm/grid/cpgrid/CpGridData.hpp>

struct PtrWrapper
{
    std::shared_ptr<Dune::cpgrid::CpGridData> data;
};
/// THis should compile
std::shared_ptr<Dune::cpgrid::CpGridData> getPointer()
{
    auto data = std::make_shared<Dune::cpgrid::CpGridData>();
    return data;
}

PtrWrapper getPointerWrapper()
{
    PtrWrapper data;
    return data;
}

std::tuple<std::shared_ptr<Dune::cpgrid::CpGridData>> getPointerTuple()
{
    auto data = std::make_shared<Dune::cpgrid::CpGridData>();
    return {data};
}

std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> getPointerVector()
{
    auto data = std::make_shared<Dune::cpgrid::CpGridData>();
    return { data };
}

Dune::cpgrid::CpGridData getCopyRVO()
{
    // This works because of copy elision
    // https://en.cppreference.com/w/cpp/language/copy_elision
    // SO no copy constructor is called
    return {};
}

struct Wrapper
{
    Dune::cpgrid::CpGridData data;
};

/// These won't compile, because copy constuctor is private/deleted

// error: ‘Dune::cpgrid::CpGridData::CpGridData(const Dune::cpgrid::CpGridData&)’ is private within this context
Dune::cpgrid::CpGridData getCopy()
{
    Dune::cpgrid::CpGridData data;
    return data;
}

// error: use of deleted function ‘Wrapper::Wrapper(const Wrapper&)’
Wrapper getWrapper()
{
    Wrapper data;
    return data;
}

// error: could not convert ‘{data}’ from ‘<brace-enclosed initializer list>’ to ‘std::tuple<Dune::
std::tuple<Dune::cpgrid::CpGridData> getTuple()
{
    Dune::cpgrid::CpGridData data;
    return {data};
}

// error: static assertion failed: result type must be constructible from value type of input range
std::vector<Dune::cpgrid::CpGridData> getVector()
{
    Dune::cpgrid::CpGridData data;
    return {data};
}

int main()
{
    auto ret1 = getPointer();
    auto ret2 = getPointerWrapper();
    auto ret3 = getPointerTuple();
    auto ret4 = getPointerVector();
    auto ret5 = getCopyRVO();
    auto ret6 = getCopy();
    auto ret7 = getWrapper();
    auto ret8 = getTuple();
    auto ret9 = getVector();
}
