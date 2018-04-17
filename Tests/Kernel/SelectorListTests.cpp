// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AbortList Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Fatras/Kernel/SelectorList.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

namespace Test {

  /// needed are : cache, generator, detector, particle
  struct Object
  {
    int         feature = 0;
    std::string name    = "";
  };  
  
  /// Selector that selectos on the Feature
  struct FeatureSelector {

    int select_on = 0;
    
    /// call operator 
    template <typename particle_t>
    bool
    operator()(const particle_t& object) const 
    { return object.feature == select_on; }
  
  };
  
  /// Selector that selectos on the Name
  struct NameSelector {

    std::string select_on = "";
    
    /// call operator 
    template <typename particle_t>
    bool
    operator()(const particle_t& object) const 
    { return object.name == select_on; }
  
  };
  
  // This tests the implementation of the selector list 
  BOOST_AUTO_TEST_CASE(SelectorList_test)
  {

    // An object with name and features
    Object o1;
    o1.name    = "o";
    o1.feature = 1;
    
    // Selector that trifers on the feature
    FeatureSelector selector1;
    selector1.select_on =1;
    // the wrong feature value
    FeatureSelector selector2;
    selector2.select_on =2;
    
    // test that the feature value 1 is selected
    BOOST_TEST(selector1(o1));
    // test that the feature value 2 is selected
    BOOST_TEST(!selector2(o1));
    
    // Let's test this with the selector list
    SelectorList<FeatureSelector> selectorList11;
    auto& sl11 = selectorList11.template get<FeatureSelector>();
    sl11.select_on = 1;
    
    SelectorList<FeatureSelector> selectorList12;
    auto& sl12 = selectorList12.template get<FeatureSelector>();
    sl12.select_on = 2;

    // test that the feature value 1 is selected
    BOOST_TEST(selectorList11(o1));
    // test that the feature value 2 is selected
    BOOST_TEST(!selectorList12(o1));
    
    // make a combined selector lsit
    SelectorList<FeatureSelector,NameSelector> o1List;
    auto& s1 = o1List.get<FeatureSelector>();
    s1.select_on = 1;
    auto& so = o1List. get<NameSelector>();
    so.select_on = "o";
    
    // test that the feature value 1 is selected
    BOOST_TEST(o1List(o1));

    // make a combined selector lsit
    SelectorList<FeatureSelector,NameSelector> o2List;
    auto& s2 = o2List.template get<FeatureSelector>();
    s2.select_on = 2;
    so = o2List.template get<NameSelector>();
    so.select_on = "o";
    
    // test that the feature value 1 is selected
    BOOST_TEST(!o2List(o1));
  }

}  // namespace Test
}  // namespace Fatras
