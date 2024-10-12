#include "BoardGameBotAI.hpp"


// Standard lib
#include <cassert>

// Algo headers

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


BoardGameBotAI::BoardState *BoardGameBotAI::CreateBoard(const Field::Field2<char> &iPawns,
                                                        const std::vector<std::array<int, 2>> &iMove,
                                                        const int iDepth) {
  BoardState *newBoard= new BoardState;
  newBoard->Pawns= iPawns;
  newBoard->Move= iMove;
  newBoard->Score= 0;
  newBoard->NashScore= (IsRedTurn(iDepth)) ? -INT_MAX : INT_MAX;
  newBoard->NashNbSteps= INT_MAX;
  return newBoard;
}


void BoardGameBotAI::DeleteBoard(BoardState *ioBoard) {
  assert(ioBoard != nullptr);
  for (BoardState *subBoard : ioBoard->SubBoards)
    DeleteBoard(subBoard);
  delete ioBoard;
  ioBoard= nullptr;
}
