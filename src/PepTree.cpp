#include "PepTree.hpp"

#include <cstdio>
#include <cassert>
#include <queue>
#include <cstdint>
#include <limits>
#include <boost/range/algorithm/for_each.hpp>

#include "Fasta.hpp"

using namespace std;
using boost::range::for_each;

size_t Trie::GetLeafCreatePath( char const * seq ) {
		size_t curIndex = 0, i = 0;

		while ( i != Depth() ) {
				auto node = GetNode( curIndex );
				// Do not use map::operator[] It creates the nodes, effectively creating the complete tree
				auto it = node->children.find( seq[i] );
				if ( it != node->children.end() ) {
						curIndex = it->second;
						++i;
				} else {
						break;
				}
		}

		if ( i != Depth() ) {
				for ( ; i != Depth() - 1; ++i ) {
						auto index = curIndex;
						curIndex = AllocateNode();   // we want curIndex to change for next iteration
						GetNode( index )->children[seq[i]] = curIndex;
				}
				auto index = curIndex;
				curIndex = AllocateLeaf();   // we want curIndex to point to the leaf index
				GetNode( index )->children[seq[i]] = curIndex;
		}

		return curIndex;
}

namespace {

	size_t LinearizeLeafData( vector< LeafBaseDataType > & arrayLeafData, Trie const & trie, size_t leafIndex ) {
			auto leaf = trie.GetLeaf( leafIndex );
			auto curIndex = arrayLeafData.size();

			if ( leaf->positions.size() > numeric_limits< uint16_t >::max() ) {
					throw std::runtime_error{ "Leaf data vector-of-proteins size overflow, abording" };
			}
			auto vecSz = static_cast< uint16_t >( leaf->positions.size() );
			arrayLeafData.push_back( static_cast< LeafBaseDataType >( vecSz >> 8 ) );
			arrayLeafData.push_back( static_cast< LeafBaseDataType >( vecSz      ) );

			for_each( leaf->positions, [&]( pair< uint32_t, vector< size_t > > const & p ) {
					auto protIdx = p.first;
					arrayLeafData.push_back( static_cast< char >( protIdx >> 24 ) );
					arrayLeafData.push_back( static_cast< char >( protIdx >> 16 ) );
					arrayLeafData.push_back( static_cast< char >( protIdx >>  8 ) );
					arrayLeafData.push_back( static_cast< char >( protIdx       ) );

					if ( p.second.size() > numeric_limits< uint16_t >::max() ) {
							throw std::runtime_error{ "Leaf data vector-of-positions size overflow, abording" };
					}
					uint16_t sz = (uint16_t)p.second.size();
					arrayLeafData.push_back( (LeafBaseDataType)(sz >> 8) );
					arrayLeafData.push_back( (LeafBaseDataType)sz );

					for_each( p.second, [&]( size_t z ) {
							if ( z > numeric_limits< uint16_t >::max() ) {
									throw std::runtime_error{ "Leaf data position overflow, abording" };
							}
							uint16_t pos = (uint16_t)z;
							arrayLeafData.push_back( (LeafBaseDataType)(pos >> 8) );
							arrayLeafData.push_back( (LeafBaseDataType)pos );
					});
			});

			return curIndex;
	}

	size_t LinearizeLeaves( vector< LeafBaseDataType > & arrayLeaves
	                      , vector< LeafBaseDataType > & arrayLeafData
	                      , Trie const & trie, size_t leafIndex
	                      , string const & str
	                      ) {
			size_t curIndex = arrayLeaves.size();
			size_t strSzIn32b = StringBufferSizeInWords( str.size() );
			arrayLeaves.resize( curIndex + (strSzIn32b + 1)*sizeof( uint32_t ) );   // reserve size for buffer + link
			std::copy( begin( str ), end( str ), arrayLeaves.data() + curIndex );
			size_t link = LinearizeLeafData( arrayLeafData, trie, leafIndex );
			if ( link > numeric_limits< uint32_t >::max() ) {
					throw std::runtime_error{ "Leaf data index overflow, abording" };
			}
			arrayLeaves[arrayLeaves.size()-4] = (LeafBaseDataType)(link >> 24);
			arrayLeaves[arrayLeaves.size()-3] = (LeafBaseDataType)(link >> 16);
			arrayLeaves[arrayLeaves.size()-2] = (LeafBaseDataType)(link >>  8);
			arrayLeaves[arrayLeaves.size()-1] = (LeafBaseDataType)link;
			return curIndex;
	}

	struct NodeQueueData {
			size_t                 nodeIndex;
			string                 genealogy;
			uint32_t             * updateParent;   // !null if first child and parent 'beginning of children' link is to be updated
			vector< uint32_t * >   leafUpdates;
	};
	static size_t const endOfChildrenSentinelIndex = ((size_t)-1);
	static size_t const rootSentinelIndex          = ((size_t)-2);

	// Complex linearization by breadth first traversal in order for the range
	// described by node n and n+1 take into account every child leave of the subtree
	void LinearizeNodes( vector< EncodedNodeType >  & arrayNodes
	                   , vector< LeafBaseDataType > & arrayLeaves
	                   , vector< LeafBaseDataType > & arrayLeafData
	                   , Trie const & trie
	                   ) {
			queue< NodeQueueData > nodeQueue;
			nodeQueue.push( NodeQueueData{ rootSentinelIndex, {}, nullptr, {} } );

			while ( !nodeQueue.empty() ) {
					auto nodeData = move( nodeQueue.front() );
					nodeQueue.pop();

					size_t thisIndex = arrayNodes.size();
					if ( nodeData.nodeIndex == endOfChildrenSentinelIndex) {   // dummy node
							arrayNodes.resize( thisIndex + 1 );   // default init to 0 = end of children sentinel
							continue;
					} else if ( nodeData.nodeIndex != rootSentinelIndex ) {
							arrayNodes.resize( thisIndex + 2 );   // size for encoded letter+childIndex AND index to first leaf

							// the order if(){...} then affectation here is important
							// since the root node have the "same index" as its first child
							// if written the other way around we end up
							// overwritting the stored aa character with an invalid compressed components parts
							if ( nodeData.updateParent ) {
									if ( thisIndex > 134217727 ) {   //< 2^27-1 ~ UINT27_MAX
											throw std::runtime_error{ "Node index overflow, abording" };
									}
									*nodeData.updateParent = CompressComponents( { (char)*nodeData.updateParent, (uint32_t)(thisIndex) } );
							}
							arrayNodes[thisIndex] = nodeData.genealogy.back();
					}

					if ( nodeData.genealogy.size() < trie.Depth() ) {
							auto node = trie.GetNode( nodeData.nodeIndex == rootSentinelIndex ? 0 : nodeData.nodeIndex );
							bool firstChild = true;
							for_each( node->children, [&]( pair< char, size_t > const & p ) {
									if ( firstChild ) {
											nodeData.leafUpdates.push_back( &arrayNodes[thisIndex + 1] );
											nodeQueue.push( NodeQueueData{ p.second
											                             , nodeData.genealogy + p.first
											                             , &arrayNodes[thisIndex]
											                             , move( nodeData.leafUpdates )
											                             } );
											firstChild = false;
									} else {
											nodeQueue.push( NodeQueueData{ p.second
											                             , nodeData.genealogy + p.first
											                             , nullptr
											                             , vector< uint32_t * >{}
											                             } );
									}
							});
							nodeQueue.push( NodeQueueData{ endOfChildrenSentinelIndex, {}, nullptr, {} } );
					} else {
							uint32_t childLeafIndex = LinearizeLeaves( arrayLeaves, arrayLeafData
							                                         , trie, nodeData.nodeIndex
							                                         , nodeData.genealogy
							                                         );
							if ( childLeafIndex > 134217727 ) {   //< 2^27-1 ~ UINT27_MAX
									throw std::runtime_error{ "Leaf index overflow, abording" };
							}

							for_each( nodeData.leafUpdates, [childLeafIndex]( uint32_t * linkPtr ) {
									*linkPtr = EncodeLeafLink( childLeafIndex );
							});
							arrayNodes[thisIndex]   = CompressComponents( { nodeData.genealogy.back(), 0 } );
							arrayNodes[thisIndex+1] = EncodeLeafLink( childLeafIndex );
					}
			}
	}

}

LinearizedTree Trie::LinearizeTree() const {
		vector< EncodedNodeType >  arrayNodes;
		vector< LeafBaseDataType > arrayLeaves;
		vector< LeafBaseDataType > arrayLeafData;

		// nb link = nb leaves + nb internal*2 (includes leaves link per internal, excluding root) = Leaf::nbLeaves + 2*(Node::nbNodes - 1),
		//   plus trailing 0s to indicate 'end of children list', which is one per node (including root) = Node::nbNodes
		arrayNodes.reserve( 2*(NumLeaves() + NumNodes() - 1) + NumNodes() );
		arrayLeaves.reserve( NumLeaves()*(StringBufferSizeInWords( Depth() ) + 1/*link*/)*sizeof( uint32_t ) + 1/*sentinel*/ );

		LinearizeNodes( arrayNodes, arrayLeaves, arrayLeafData, *this );

		arrayLeaves.push_back( '\0' );   // end of leaves sentinel

		return LinearizedTree{ Depth(), arrayNodes, arrayLeaves, arrayLeafData };
}

static uint32_t const nodesOffset = 3 * sizeof( uint32_t );   // treeDepth + LeavesOffset + LeafDataOffset

void WriteLinearizedTree( FILE * file, LinearizedTreeData const & td ) {
		uint32_t     treeDepth     = GetTreeDepth( td );

		auto const & nodesVector   = GetTreeNodes( td );
		auto const & leavesVector  = GetTreeLeaves( td );
		auto const & leafPosVector = GetTreeLeafPos( td );

		size_t nodesSize   = GetVectorSize( nodesVector );
		size_t leavesSize  = GetVectorSize( leavesVector );
		size_t leafPosSize = GetVectorSize( leafPosVector );

		EncodedNodeType const  * nodesData   = GetVectorData( nodesVector );
		LeafBaseDataType const * leavesData  = GetVectorData( leavesVector );
		LeafBaseDataType const * leafPosData = GetVectorData( leafPosVector );

		uint32_t leavesOffset  = nodesOffset  + nodesSize  * sizeof( *nodesData );
		uint32_t leafPosOffset = leavesOffset + leavesSize * sizeof( *leavesData );

		uint32_t arr[] = { treeDepth, leavesOffset, leafPosOffset };
		fwrite( arr        , sizeof( arr[0] )      , sizeof( arr ) / sizeof( arr[0] ), file );
		fwrite( nodesData  , sizeof( *nodesData )  , nodesSize                       , file );
		fwrite( leavesData , sizeof( *leavesData ) , leavesSize                      , file );
		fwrite( leafPosData, sizeof( *leafPosData ), leafPosSize                     , file );
}

void WriteReadableLinearizedNodes( FILE * file, LinearizedTreeData const & td ) {
		uint32_t index = 0;
		size_t nodeNumber = 0;
		uint32_t const * p   = GetTreeNodesData( td );
		uint32_t const * end = p + GetTreeNodesSize( td );
		while( p != end ) {
				uint32_t val = *p++;
				auto comp = ExtractComponents( val );
				char c = GetAAChar( comp );
				if ( IsEncodedLeafLink( val ) ) {
						fprintf( file, "        %06X: --> %06X\n", index++, ExtractEncodedLeafLink( val ) );
				} else if ( Fasta::IsValidAA( c ) ) {
						fprintf( file, "(%05zu) %06X: %c %06X\n", nodeNumber++, index++, c, GetChildIndex( comp ) );
				} else {
						fprintf( file, "        %06X: |\n", index++ );
				}
		}
}

void WriteReadableLinearizedLeaves( FILE * file, LinearizedTreeData const & td ) {
		uint32_t index = 0;
		size_t nodeNumber = 0;
		ForEachLeaf( td, [&]( char const * str, uint32_t dataOffset ) {
				printf( "(%05zu) %06X: %s -> %06X\n", nodeNumber++, index, str, dataOffset );
				index += (StringBufferSizeInWords( GetTreeDepth( td ) ) + 1)*sizeof( uint32_t );
		});
}

namespace {

	struct ReadablePrinterFunctor {
			bool firstPos;
			FILE * file;
			ReadablePrinterFunctor( FILE * f ): firstPos( true ), file( f ) { }

			void ListSize( uint16_t ) {   }

			void AddHeader( uint16_t protIndex, uint16_t s ) {
					fprintf( file, " %06hX[", protIndex );
			}
			void StopHeader() {   }

			void AddPos( uint16_t p ) {
					if ( !firstPos ) {
							fprintf( file, ", " );
					}
					fprintf( file, "%hu", p );
					firstPos = false;
			}
			void StopPos() {
					fprintf( file, "]" );
					firstPos = true;
			}
	};

}

void WriteReadableLinearizedLeafPos( FILE * file, LinearizedTreeData const & td ) {
		uint32_t index = 0;
		size_t nodeNumber = 0;
		while ( index < GetTreeLeafPosSize( td ) ) {
				fprintf( file, "(%05zu) %06tX:", nodeNumber++, index );
				index += ForLeafPos( td, index, ReadablePrinterFunctor{ file } );
				fprintf( file, "\n" );
		}
}

void WriteReadableLinearizedTree( FILE * file, LinearizedTreeData const & td ) {
		WriteReadableLinearizedNodes( file, td );
		fprintf( file, "~~~~~~\n" );
		WriteReadableLinearizedLeaves( file, td );
		fprintf( file, "~~~~~~\n" );
		WriteReadableLinearizedLeafPos( file, td );
}

LinearizedTree ReadLinearizedTree( FILE * file ) {
		fseek( file, 0, SEEK_END );
		auto fileSize = ftell( file );
		fseek( file, 0, SEEK_SET );

		uint32_t arr[3];
		fread( arr, sizeof( arr[0] ), sizeof( arr ) / sizeof( arr[0] ), file );

		vector< EncodedNodeType > nodes( (arr[1] - nodesOffset)/sizeof( EncodedNodeType ) );
		fread(  nodes.data(),    sizeof( nodes[0] ), nodes.size(), file );
		
		vector< LeafBaseDataType > leaves( arr[2] - arr[1] );
		fread( leaves.data(),   sizeof( leaves[0] ), leaves.size(), file );

		vector< LeafBaseDataType > leafData( fileSize - arr[2] );
		fread( leafData.data(), sizeof( leafData[0] ), leafData.size(), file );

		return LinearizedTree{ arr[0], nodes, leaves, leafData };
}

size_t GetTreeNumberLeaves( FILE * file ) {
		fseek( file, 0, SEEK_END );
		ftell( file );
		fseek( file, 0, SEEK_SET );

		uint32_t arr[3];
		fread( arr, sizeof( arr[0] ), sizeof( arr ) / sizeof( arr[0] ), file );

		return (arr[2] - arr[1]) / LeavesLinkSize( arr[0] );
}

LinearizedTreeData ResolveLinearizedTree( Byte const * bytes, size_t size ) {
		uint32_t depth         = (bytes[0] << 24) + (bytes[1] << 16) + (bytes[ 2] << 8) + bytes[ 3];
		uint32_t leavesOffset  = (bytes[3] << 24) + (bytes[5] << 16) + (bytes[ 6] << 8) + bytes[ 7];
		uint32_t leafPosOffset = (bytes[8] << 24) + (bytes[9] << 16) + (bytes[10] << 8) + bytes[11];
		return LinearizedTreeData{ depth
		                         , { (EncodedNodeType  *)(bytes + nodesOffset)  , leavesOffset  - nodesOffset   }
		                         , { (LeafBaseDataType *)(bytes + leavesOffset) , leafPosOffset - leavesOffset  }
		                         , { (LeafBaseDataType *)(bytes + leafPosOffset), size          - leafPosOffset }
		                         };
}
