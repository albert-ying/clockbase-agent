import os
import uvicorn
from dotenv import load_dotenv
from fastapi.middleware.cors import CORSMiddleware
from fastapi import FastAPI, Header, HTTPException, Depends, status
from pydantic import BaseModel, field_validator
from typing import Optional


from paper_qa.pipe import PaperQAPipe

load_dotenv()

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# A simple function to simulate token verification (replace with your actual logic)
def verify_token(token: str) -> bool:
    # In a real application, you would:
    # 1. Decode the token (e.g., using JWT libraries like `python-jose`)
    # 2. Validate its signature
    # 3. Check its expiry
    # 4. Potentially retrieve user information from the token payload
    print(f"Verifying token: {token}")
    return token == os.getenv("LOCAL_TOKEN") # For demonstration purposes only!

# Dependency to extract and validate the bearer token
async def get_bearer_token(authorization: Optional[str] = Header(None)) -> str:
    if authorization is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Authorization header is missing",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    # Check if the header starts with "Bearer "
    if not authorization.startswith("Bearer "):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid Authorization header format. Must be 'Bearer <token>'",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    token = authorization.split(" ")[1] # Extract the token part
    
    if not verify_token(token):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or expired token",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    return token

class QARequest(BaseModel):
    token: str
    user_input: str

    @field_validator('token')
    def validate_token(cls, v):
        if not verify_token(v):
            raise ValueError("Invalid or expired token")
        return v

class QAResponse(BaseModel):
    message: str
    status: str = 'success'

@app.post("/qa")
async def qa(request: QARequest):
    pipe = PaperQAPipe()
    try:
        answer = await pipe.pipeline(request.user_input)
        return QAResponse(message=answer, status='success')
    except Exception as e:
        return QAResponse(message=str(e), status='error')


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.options("/qa")
async def options():
    return {"status": "ok"}


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)