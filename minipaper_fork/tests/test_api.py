import requests
import os
from dotenv import load_dotenv

load_dotenv()


def test_api():
    url = "http://localhost:8000/qa"
    headers = {
        "Authorization": f"Bearer {os.getenv('LOCAL_TOKEN')}"
    }
    data = {"user_input": "Frequently mentioned potential “routine laboratory” biomarkers"}
    response = requests.post(url, json=data, headers=headers)
    print(response.json())


if __name__ == "__main__":
    test_api()