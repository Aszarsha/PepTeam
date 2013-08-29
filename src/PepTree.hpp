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
#define BITS_FOR_LINK_MASK ((unsigned)(1<<(NB_BITS_FOR_LINK))-1)
#define BITS_FOR_CHAR_MASK ~BITS_FOR_LINK_MASK

typedef unsigned char Byte;
typedef uint32_t      EncodedNodeType;
typedef Byte          LeafBaseDataType;

inline uint32_t StringBufferSizeInWords( size_t strSz ) {
		size_t bufSz = strSz + 1;   // For \0
		return bufSz / sizeof( uint32_t ) + ((bufSz % sizeof( uint32_t )) != 0 ? 1 : 0);
}

inline uint32_t LeavesLinkSize( size_t treeDepth ) {
		return (StringBufferSizeInWords( treeDepth ) + 1)*sizeof( uint32_t );
}

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

class MemPepTree {
	public:
		MemPepTree( uint32_t depth
		          , std::vector< EncodedNodeType >  && nodes
		          , std::vector< LeafBaseDataType > && leaves
		          , std::vector< LeafBaseDataType > && leafPos
		          );

	public:
		void Write( FILE * file ) const;

	private:
		uint32_t depth;
		std::vector< EncodedNodeType >  nodes;
		std::vector< LeafBaseDataType > leaves;
		std::vector< LeafBaseDataType > leafPos;
};

class MMappedPepTree {
	public:
		explicit MMappedPepTree( const char * filename );

		~MMappedPepTree();

	public:
		uint32_t Depth() const;

		uint32_t GetNumberLeaves() const;

		void WriteReadableTree   ( FILE * file ) const;
		void WriteReadableNodes  ( FILE * file ) const;
		void WriteReadableLeaves ( FILE * file ) const;
		void WriteReadableLeafPos( FILE * file ) const;

		EncodedNodeType const * GetNodesData() const;
		size_t                  GetNodesSize() const;

		LeafBaseDataType const * GetLeavesData() const;
		size_t                   GetLeavesSize() const;

		LeafBaseDataType const * GetLeafPosData() const;
		size_t                   GetLeafPosSize() const;

	public:
		template< typename F >
		inline void ForNodeChildren( uint32_t index, F && f ) const {
				auto data = GetNodesData();
				size_t childNumber = 0;
				while ( true ) {
						EncodedNodeType val = data[index++];
						if ( val == 0 ) {
								break;
						}

						EncodedNodeType leafStart = ExtractEncodedLeafLink( data[index++] );
						EncodedNodeType leafStop;
						if ( index < GetNodesSize() - 1 ) {
								leafStop = data[index+1];
								if ( data[index] == 0 ) {   // presence of end-of-children sentinel
										leafStop = data[index+2];
								}
								leafStop = ExtractEncodedLeafLink( leafStop  );
								if ( leafStop < leafStart ) {   // new tree depth line
										leafStop = GetLeavesSize();
								}
						} else {
								leafStop = GetLeavesSize();
						}
						auto comp = ExtractComponents( val );
						f( childNumber, GetAAChar( comp ), GetChildIndex( comp ), leafStart, leafStop );
						++childNumber;
				}
		}

		template< typename F >
		inline void ForLeaf( uint32_t index, F && f ) const {
				auto data = GetLeavesData() + (index * LeavesLinkSize( Depth() ));
				size_t alignedBufSz = StringBufferSizeInWords( Depth() );
				ExtractLeafCallF( data, alignedBufSz, std::forward< F >( f ) );
		}

		template< typename F >
		inline void ForLeafRange( uint32_t start, uint32_t stop, F && f ) {
				auto data = GetLeavesData() + (start * LeavesLinkSize( Depth() ));
				size_t alignedBufSz = StringBufferSizeInWords( Depth() );
				while ( start < stop ) {
						ExtractLeafCallF( data, alignedBufSz, std::forward< F >( f ) );
						start += alignedBufSz + 1;
				}
		}

		template< typename F >
		inline void ForEachLeaf( F && f ) const {
				auto data = GetLeavesData();
				size_t alignedBufSz = StringBufferSizeInWords( Depth() );
				while ( *data != '\0' ) {
						ExtractLeafCallF( data, alignedBufSz, std::forward< F >( f ) );
						data += (alignedBufSz + 1)*sizeof( uint32_t );
				}
		}

		template< typename F >
		inline uint32_t ForLeafPos( uint32_t index, F && f ) const {
				auto base = GetLeafPosData() + index;
				auto data = base;

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

	private:
		int fd;
		uint32_t fileSize;
		char const * ptr;

	private:
		template< typename F >
		inline void ExtractLeafCallF( Byte const * data, size_t alignedBufSz, F && f ) const {
				Byte const * p = data + alignedBufSz*sizeof( uint32_t );
				uint32_t offset = *p++ << 24;
				offset         += *p++ << 16;
				offset         += *p++ <<  8;
				offset         += *p;
				f( (char const *)data, offset );
		}
};

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
		explicit Trie( uint32_t treeDepth )
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

		uint32_t Depth() const {   return depth;   }

		size_t NumNodes() const  {   return nodes .size();    }
		size_t NumLeaves() const {   return leaves.size();    }

		size_t GetLeafCreatePath( char const * seq );

		MemPepTree LinearizeTree() const;

	private:
		uint32_t depth;

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
