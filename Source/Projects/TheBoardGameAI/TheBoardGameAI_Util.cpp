#include "TheBoardGameAI.hpp"


// Standard lib

// Algo headers

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


TheBoardGameAI::BoardState *TheBoardGameAI::CreateBoard(const Field::Field2<int> &iPawns,
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


void TheBoardGameAI::DeleteBoard(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] DeleteBoard on a null board\n");
  for (BoardState *subBoard : ioBoard->SubBoards)
    DeleteBoard(subBoard);
  delete ioBoard;
  ioBoard= nullptr;
}
