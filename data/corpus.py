import json
import os
import glob

if __name__ == "__main__":

    print("running")
    # Parse Pyrrha data to files
    files = glob.glob("lemmat_src/*")

    data = []

    for file in files:
        print("processing " + file)

        with open(file, 'r') as f:
            f.readline()
            myText = {"textID": file.split("/")[1].rstrip().split("-")[0],
                      "content": []
                      }
            for line in f.readlines():
                #TODO: correct that
                if len(line.split("\t")) >= 3:
                    myText["content"].append(
                        {"token": line.split('\t')[0], "lemma": line.split('\t')[1], "pos": line.split('\t')[2]}
                    )

        data.append(myText)

    for doc in data:
        with open("LEMMA/"+doc["textID"], 'w') as out:
            for line in doc["content"]:
                out.write(line["lemma"] + " ")

    for doc in data:
        with open("POS/"+doc["textID"], 'w') as out:
            for line in doc["content"]:
                if not line["pos"].startswith("PON"):
                    out.write(line["pos"] + " ")


    for doc in data:
        with open("TXT/"+doc["textID"], 'w') as out:
            for line in doc["content"]:
                out.write(line["token"] + " ")


    myjson = open("function_words_old-french_of3c.json", "r")
    fw = [f[0] for f in json.load(myjson)]

    for doc in data:
        with open("MOTIFS/"+doc["textID"], 'w') as out:
            for line in doc["content"]:
                if line["lemma"] in fw:
                    out.write(line["lemma"] + " ")
                else:
                    if not line["pos"].startswith("PON"):
                        out.write(line["pos"] + " ")


