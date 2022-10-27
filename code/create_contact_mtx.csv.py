import csv

row_list = []
for row in range(15):
    new_row = []
    for col in range(15):
        if row == 2 and col == 2:
            new_row.append(1)
        elif row == 4 and col == 4:
            new_row.append(1)
        # elif row == 6 and col == 6:
        #     new_row.append(1)
        # elif row == 8 and col == 8:
        #     new_row.append(1)
        else:
            new_row.append(0)
    row_list.append(new_row)


with open('contact_mtx_2_small.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(row_list)

