#ifndef DETECTOR_DESCRIPTION_CORE_DD_COMPACT_VIEW_H
#define DETECTOR_DESCRIPTION_CORE_DD_COMPACT_VIEW_H

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "DetectorDescription/Core/interface/DDCompactViewImpl.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"
#include "DetectorDescription/Core/interface/DDTranslation.h"
#include "DetectorDescription/Core/interface/Store.h"
#include "DetectorDescription/Core/interface/DDLogicalPart.h"
#include "DetectorDescription/Core/interface/DDPosData.h"
#include "DetectorDescription/Core/interface/DDTransform.h"
#include "DataFormats/Math/interface/Graph.h"
#include "DataFormats/Math/interface/GraphWalker.h"

class DDDivision;
class DDName;
struct DDPosData;

/**
  Navigation through the compact view of the detector ...

Updated: Michael Case [ MEC ] 2010-02-11

*/
//MEC: these comments are kept from original... Will we ever do this? don't think so.
//FIXME: DDCompactView: check for proper acyclic directed graph structure!!
//FIXME:
//FIXME:         X          [A-X] ... LogicalPart
//FIXME:        / \             | ... PosPart (directed parten to child)
//FIXME:       A   A
//FIXME:       |   | 
//FIXME:       B   C      
//FIXME:
//FIXME:    THIS IS NOT ALLOWED, but currently can be specified using DDL ....
//FIXME:

//! Compact representation of the geometrical detector hierarchy
/** A DDCompactView represents the detector as an acyclic directed multigraph.
    The nodes are instances of DDLogicalPart while the edges are pointers to
    DDPosData. Edges are directed from parent-node to child-node. 
    Edges represented by DDPosData are the relative translation and rotation
    accompanied by a copy-number of the child-node towards the parent-node.
    
    One node is explicitly marked as the root node. It is the DDLogicalPart which
    defines the base coordinate system for the detector. All subdetectors are
    directly or inderectly positioned inside the root-node. 
    
    Example:
    
    The figureshows a compact-view graph consisting of 16 DDLogicalParts 
    interconnected by 20 edges represented by pointers to DDPosData.
    \image html compact-view.gif
    \image latex compact-view.eps
    
    The compact-view also serves as base for calculating nodes in the expanded
    view. Paths through the compact-view can be viewed as nodes in the expanded-view
    (expansion of an acyclic directed multigraph into a tree). In the figure there are
    5 different paths from CMS to Module2 (e.g. CMS-Pos1->Ecal-Pos4->EEndcap-Pos21->Module2)
    thus there will be 5 nodes of Module2 in the expanded view.
*/
class DDCompactView
{ 
public:
  using Graph = math::Graph<DDLogicalPart, DDPosData* >;
  using GraphWalker = math::GraphWalker<DDLogicalPart, DDPosData* >;
  
  //! Creates a compact-view 
  explicit DDCompactView();
  
  //! \b EXPERIMENTAL! Creates a compact-view using a different root of the geometrical hierarchy
  explicit DDCompactView(const DDLogicalPart & rootnodedata);
  
  //! Provides read-only access to the data structure of the compact-view.
  const Graph & graph() const;
  GraphWalker walker() const;

  //! returns the DDLogicalPart representing the root of the geometrical hierarchy
  const DDLogicalPart & root() const;
  
  //! The absolute position of the world
  const DDPosData * worldPosition() const;

  void position (const DDLogicalPart & self,
		 const DDLogicalPart & parent,
		 const std::string& copyno,
		 const DDTranslation & trans,
		 const DDRotation & rot,
		 const DDDivision * div = nullptr);
  
  void position (const DDLogicalPart & self,
		 const DDLogicalPart & parent,
		 int copyno,
		 const DDTranslation & trans,
		 const DDRotation & rot,
		 const DDDivision * div = nullptr);
  
  // ************************************************************************
  // UNSTABLE STUFF below! DON'T USE!
  // ************************************************************************
  
  //! \b don't \b use : interface not stable ....
  void setRoot(const DDLogicalPart & root);

  void lockdown();
  
 private:
  std::unique_ptr<DDCompactViewImpl> rep_;
  std::unique_ptr<DDPosData> worldpos_ ;
};

#endif
