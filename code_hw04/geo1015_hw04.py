import json, sys
from my_code_hw04 import detect_planes

def main():
    try:
        jparams = json.load(open('params.json'))
    except:
        print("ERROR: something is wrong with the params.json file.")
        sys.exit()

    detect_planes(jparams)

if __name__ == "__main__":
    main()