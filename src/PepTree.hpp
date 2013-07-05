#ifndef PEPTEAMTREE_HPP
#define PEPTEAMTREE_HPP

#include <cstdio>
#include <cstdint>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <sstream>

#include "Fasta.hpp"

// ~~~ Vector Based Tree ~~~ //
#define NB_BITS_FOR_LINK 27U
//#define BITS_FOR_LINK_MASK ((-1U)>>(32-NB_BITS_FOR_LINK))
#define BITS_FOR_LINK_MASK ((1<<(NB_BITS_FOR_LINK))-1)
#define BITS_FOR_CHAR_MASK ~BITS_FOR_LINK_MASK

typedef unsigned char Byte;
typedef uint32_t      EncodedNodeType;
typedef Byte          LeafBaseDataType;
typedef std::pair< EncodedNodeType  const *, size_t > NodesVector;
typedef std::pair< LeafBaseDataType const *, size_t > LeavesVector;
typedef std::pair< LeafBaseDataType const *, size_t > LeafPosVector;
typedef std::tuple< uint32_t
                  , std::vector< EncodedNodeType >
                  , std::vector< LeafBaseDataType >
                  , std::vector< LeafBaseDataType >
                  > LinearizedTree;
typedef std::tuple< uint32_t
                  , NodesVector
                  , LeavesVector
                  , LeafPosVector
                  > LinearizedTreeData;

template< typename V > inline auto GetVectorData( V const & v ) -> decltype( std::get< 0 >( v ) ) {   return std::get< 0 >( v );   }
template< typename V > inline auto GetVectorSize( V const & v ) -> decltype( std::get< 1 >( v ) ) {   return std::get< 1 >( v );   }

inline uint32_t                                GetTreeDepth(   LinearizedTree const     & t ) {   return std::get< 0 >( t );   }
inline uint32_t                                GetTreeDepth(   LinearizedTreeData const & t ) {   return std::get< 0 >( t );   }
inline std::vector< EncodedNodeType > const  & GetTreeNodes(   LinearizedTree const     & t ) {   return std::get< 1 >( t );   }
inline NodesVector const                     & GetTreeNodes(   LinearizedTreeData const & t ) {   return std::get< 1 >( t );   }
inline std::vector< LeafBaseDataType > const & GetTreeLeaves(  LinearizedTree const     & t ) {   return std::get< 2 >( t );   }
inline LeavesVector const                    & GetTreeLeaves(  LinearizedTreeData const & t ) {   return std::get< 2 >( t );   }
inline std::vector< LeafBaseDataType > const & GetTreeLeafPos( LinearizedTree const     & t ) {   return std::get< 3 >( t );   }
inline LeafPosVector const                   & GetTreeLeafPos( LinearizedTreeData const & t ) {   return std::get< 3 >( t );   }

inline EncodedNodeType const * GetTreeNodesData( LinearizedTree const     & t ) {   return GetTreeNodes( t ).data();             }
inline EncodedNodeType const * GetTreeNodesData( LinearizedTreeData const & t ) {   return GetVectorData( GetTreeNodes( t ) );   }
inline size_t GetTreeNodesSize( LinearizedTree const     & t ) {   return GetTreeNodes( t ).size();             }
inline size_t GetTreeNodesSize( LinearizedTreeData const & t ) {   return GetVectorSize( GetTreeNodes( t ) );   }

inline LeafBaseDataType const * GetTreeLeavesData( LinearizedTree const     & t ) {   return GetTreeLeaves( t ).data();             }
inline LeafBaseDataType const * GetTreeLeavesData( LinearizedTreeData const & t ) {   return GetVectorData( GetTreeLeaves( t ) );   }
inline size_t GetTreeLeavesSize( LinearizedTree const     & t ) {   return GetTreeLeaves( t ).size();             }
inline size_t GetTreeLeavesSize( LinearizedTreeData const & t ) {   return GetVectorSize( GetTreeLeaves( t ) );   }

inline LeafBaseDataType const * GetTreeLeafPosData( LinearizedTree const     & t ) {   return GetTreeLeafPos( t ).data();             }
inline LeafBaseDataType const * GetTreeLeafPosData( LinearizedTreeData const & t ) {   return GetVectorData( GetTreeLeafPos( t ) );   }
inline size_t GetTreeLeafPosSize( LinearizedTree const     & t ) {   return GetTreeLeafPos( t ).size();   }
inline size_t GetTreeLeafPosSize( LinearizedTreeData const & t ) {   return GetVectorSize( GetTreeLeafPos( t ) );   }

inline LinearizedTreeData GetTreeData( LinearizedTree const & t ) {
		return LinearizedTreeData{ GetTreeDepth( t )
		                         , { GetTreeNodesData( t )  , GetTreeNodesSize( t )   }
		                         , { GetTreeLeavesData( t ) , GetTreeLeavesSize( t )  }
		                         , { GetTreeLeafPosData( t ), GetTreeLeafPosSize( t ) }
		                         };
}

void WriteLinearizedTree( FILE * file, LinearizedTreeData const & t );

size_t GetTreeNumberLeaves( FILE * file );

void WriteReadableLinearizedTree( FILE * file, LinearizedTreeData const & td );
void WriteReadableLinearizedNodes( FILE * file, LinearizedTreeData const & td );
void WriteReadableLinearizedLeaves( FILE * file, LinearizedTreeData const & td );
void WriteReadableLinearizedLeafPos( FILE * file, LinearizedTreeData const & td );

LinearizedTree     ReadLinearizedTree( FILE * file );
LinearizedTreeData ResolveLinearizedTree( Byte const * bytes, size_t size );

typedef std::pair< char, uint32_t > Components;
inline char     GetAAChar( Components const & comp )     {   return std::get< 0 >( comp );   }
inline uint32_t GetChildIndex( Components const & comp ) {   return std::get< 1 >( comp );   }
inline EncodedNodeType CompressComponents( Components const & comp ) {
		uint32_t charPart = (Fasta::Char2Index( GetAAChar( comp ) )+1) << NB_BITS_FOR_LINK;
		uint32_t linkPart = GetChildIndex( comp );   // & BITS_FOR_LINK_MASK;
		return charPart | linkPart;
}
inline Components ExtractComponents( EncodedNodeType val ) {
		return { Fasta::Index2Char( (val>>NB_BITS_FOR_LINK)-1 ), val & BITS_FOR_LINK_MASK };
}
inline EncodedNodeType EncodeLeafLink( uint32_t link ) {
		uint32_t charPart = BITS_FOR_CHAR_MASK;
		uint32_t linkPart = link;   // & BITS_FOR_LINK_MASK;
		return charPart | linkPart;
}
inline bool IsEncodedLeafLink( EncodedNodeType val ) {
		return (val & BITS_FOR_CHAR_MASK) == BITS_FOR_CHAR_MASK;
}
inline uint32_t ExtractEncodedLeafLink( EncodedNodeType val ) {
		return val & BITS_FOR_LINK_MASK;
}

template< typename F >
inline void ForNodeChildren( LinearizedTreeData const & tree, uint32_t index, F && f ) {
		EncodedNodeType const * data = GetTreeNodesData( tree );
		size_t childNumber = 0;
		while ( true ) {
				EncodedNodeType val = data[index++];
				if ( val == 0 ) {
						break;
				}

				EncodedNodeType leafStart = ExtractEncodedLeafLink( data[index++] );
				EncodedNodeType leafStop;
				if ( index < GetTreeNodesSize( tree ) - 1 ) {
						leafStop = data[index+1];
						if ( data[index] == 0 ) {   // presence of end-of-children sentinel
								leafStop = data[index+2];
						}
						leafStop = ExtractEncodedLeafLink( leafStop  );
						if ( leafStop < leafStart ) {   // new tree depth line
								leafStop = GetTreeLeavesSize( tree ) - 1;   // do not count last '\0' in leaves
						}
				} else {
						leafStop = GetTreeLeavesSize( tree ) - 1;   // do not count last '\0' in leaves
				}
				auto comp = ExtractComponents( val );
				f( childNumber, GetAAChar( comp ), GetChildIndex( comp ), leafStart, leafStop );
				++childNumber;
		}
}

inline uint32_t StringBufferSizeInWords( size_t strSz ) {
		size_t bufSz = strSz + 1;   // For \0
		return bufSz / sizeof( uint32_t ) + ((bufSz % sizeof( uint32_t )) != 0 ? 1 : 0);
}

inline uint32_t LeavesLinkSize( size_t treeDepth ) {
		return (StringBufferSizeInWords( treeDepth ) + 1)*sizeof( uint32_t );
}

template< typename F >
inline void ExtractLeafCallF( Byte const * data, size_t alignedBufSz, F && f ) {
		Byte const * p = data + alignedBufSz*sizeof( uint32_t );
		uint32_t offset = *p++ << 24;
		offset         += *p++ << 16;
		offset         += *p++ <<  8;
		offset         += *p;
		f( (char const *)data, offset );
}

template< typename F >
inline void ForLeaf( LinearizedTreeData const & tree, uint32_t index, F && f ) {
		LeafBaseDataType const * data = GetTreeLeavesData( tree ) + index;
		size_t alignedBufSz = StringBufferSizeInWords( GetTreeDepth( tree ) );
		ExtractLeafCallF( data, alignedBufSz, std::forward< F >( f ) );
}

template< typename F >
inline void ForEachLeaf( LinearizedTreeData const & tree, F && f ) {
		LeafBaseDataType const * data = GetTreeLeavesData( tree );
		size_t alignedBufSz = StringBufferSizeInWords( GetTreeDepth( tree ) );
		while ( *data != '\0' ) {
				ExtractLeafCallF( data, alignedBufSz, std::forward< F >( f ) );
				data += (alignedBufSz + 1)*sizeof( uint32_t );
		}
}

template< typename F >
inline uint32_t ForLeafPos( LinearizedTreeData const & tree, uint32_t index, F && f ) {
		LeafBaseDataType const * base = GetTreeLeafPosData( tree ) + index;
		LeafBaseDataType const * data = base;

		uint16_t listSize = *data++ <<  8;
		listSize         += *data++;
		f.ListSize( listSize );
		for( uint16_t i = 0; i < listSize; ++i ) {
				uint32_t protIndex = *data++ << 24;
				protIndex         += *data++ << 16;
				protIndex         += *data++ <<  8;
				protIndex         += *data++;

				uint16_t nb = (*data++) << 8;
				nb         += *data++;
				f.AddHeader( protIndex, nb );

				for ( uint16_t i = 0; i != nb; ++i ) {
						uint16_t val = (*data++) << 8;
						val         += *data++;
						f.AddPos( val );
				}
				f.StopPos();
		}
		f.StopHeader();
		return (uint32_t)(data - base);
}

// ~~~ Node Based Tree ~~~ //
struct Trie {
	public:
		struct Node {
				std::map< char, size_t > children;
		};

		struct Leaf {
				std::map< uint32_t, std::vector< size_t > > positions;
		};

	public:
		explicit Trie( size_t treeDepth )
			: depth( treeDepth ) {
				AllocateNode();   // root
		}

	public:
		Node       * GetNode( size_t i )       {   return &nodes[i];   }
		Node const * GetNode( size_t i ) const {   return &nodes[i];   }

		Leaf       * GetLeaf( size_t i )       {   return &leaves[i];   }
		Leaf const * GetLeaf( size_t i ) const {   return &leaves[i];   }

		Node       * Root()       {   return GetNode( 0 );   }
		Node const * Root() const {   return GetNode( 0 );   }

		size_t Depth() const {   return depth;   }

		size_t NumNodes() const  {   return nodes .size();    }
		size_t NumLeaves() const {   return leaves.size();    }

		size_t GetLeafCreatePath( char const * seq );

		LinearizedTree LinearizeTree() const;

	private:
		size_t depth;

		std::vector< Node > nodes;
		std::vector< Leaf > leaves;

	private:
		size_t AllocateNode() {
				nodes.emplace_back();
				return nodes.size()-1;
		}

		size_t AllocateLeaf() {
				leaves.emplace_back();
				return leaves.size()-1;
		}
};

#endif
