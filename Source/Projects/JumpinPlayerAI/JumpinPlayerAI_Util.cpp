#include "JumpinPlayerAI.hpp"


// Standard lib

// Algo headers

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


JumpinPlayerAI::BoardState *JumpinPlayerAI::CreateBoard(const std::vector<std::vector<int>> &iPawns,
                                                        const std::array<int, 4> &iMove,
                                                        const int iDepth) {
  BoardState *newBoard= new BoardState;
  newBoard->Pawns= iPawns;
  newBoard->Move= iMove;
  newBoard->Score= (IsRedTurn(iDepth)) ? -INT_MAX : INT_MAX;
  newBoard->NashScore= newBoard->Score;
  newBoard->NashNbSteps= INT_MAX;
  return newBoard;
}


void JumpinPlayerAI::DeleteBoard(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] DeleteBoard on a null board\n");
  for (BoardState *subBoard : ioBoard->SubBoards)
    DeleteBoard(subBoard);
  delete ioBoard;
  ioBoard= nullptr;
}
